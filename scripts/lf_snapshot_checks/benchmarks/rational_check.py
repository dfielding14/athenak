"""Small CPU checks of positive rational LF filters in a curved periodic field."""
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy.sparse.linalg import LinearOperator, cg

# Fit the continuum inverse absolute derivative over a stated band.
kfit = np.geomspace(.25, 32, 401)
kcheck = np.geomspace(.25, 32, 10001)
poles = np.geomspace(.125, 64, 8)
response = kfit[:, None] * poles / (kfit[:, None]**2 + poles**2)
one = np.ones(kfit.size)
fit = linprog(np.r_[np.zeros(8), 1.],
    A_ub=np.vstack((np.column_stack((response, -one)),
                    np.column_stack((-response, -one)))),
    b_ub=np.r_[one, -one], bounds=(0, None), method='highs')
assert fit.success
weights = fit.x[:-1] * poles
fiterror = np.max(np.abs(np.sum(kcheck[:,None]*weights /
    (kcheck[:,None]**2+poles**2), axis=1)-1))
results = {'band': [.25,32], 'poles': poles.tolist(), 'weights': weights.tolist(),
           'continuum_max_relative_error': float(fiterror), 'cases': []}
for n in (16, 24):
    # Odd Fourier derivative: explicitly zero even-grid Nyquist. Its nullspace
    # includes Nyquist modes, so the benchmark uses smooth resolved modes only.
    coords = np.arange(n)*2*np.pi/n
    z,y,x = np.meshgrid(coords,coords,coords,indexing='ij')
    a = 1.
    bx, by, bz = a*np.cos(z)/np.sqrt(2), a*np.sin(z)/np.sqrt(2), 1/np.sqrt(2)
    wave = np.fft.fftfreq(n)*n
    wave[n//2] = 0
    def deriv(f, axis):
        shape = [1,1,1]; shape[axis] = n
        return np.fft.ifft(1j*wave.reshape(shape)*np.fft.fft(f,axis=axis),axis=axis).real
    def direction(f):
        return bx*deriv(f,2)+by*deriv(f,1)+bz*deriv(f,0)
    def positive(v):
        return -direction(direction(v.reshape((n,n,n)))).ravel()
    for mode in (0,1,3, 'mixed'):
        phase = x-a*np.sin(z)+(mode if mode != 'mixed' else 1)*z
        temp = np.cos(phase)
        if mode == 'mixed':
            temp = np.random.default_rng(31415).normal(size=temp.shape)
        driver = direction(temp)
        kpar = mode/np.sqrt(2) if mode != 'mixed' else None
        exact = np.sin(phase) if mode else np.zeros_like(temp)
        flux = np.zeros_like(temp)
        iterations = []
        start=time.perf_counter()
        for pole,weight in zip(poles,weights):
            operator=LinearOperator((n**3,n**3),matvec=lambda v,p=pole: positive(v)+p*p*v)
            count=[0]
            def tick(v): count[0]+=1
            sol,status=cg(operator,driver.ravel(),rtol=2e-8,atol=1e-12,maxiter=1500,callback=tick)
            assert status==0, (n,mode,pole,status)
            flux -= weight*sol.reshape(temp.shape)
            iterations.append(count[0])
        rms=lambda v: float(np.sqrt(np.mean(v*v)))
        error=(rms(flux-exact)/(rms(exact) if mode else 1.)) if mode != 'mixed' else None
        dissipation = float(-np.mean(flux*driver))
        assert dissipation >= -1e-14
        if mode in (1,3):
            assert error < .002
        results['cases'].append({'n':n,'mode':mode,'kparallel':kpar,
            'minus_mean_flux_times_gradient': dissipation,
            'rms_error_relative_or_null_absolute':error,
            'driver_rms':rms(driver),'cg_iterations_by_pole':iterations,
            'cpu_seconds':time.perf_counter()-start})
        print(results['cases'][-1],flush=True)
results['limitations'] = ['Scalar constant-density, constant-strength field only.',
 'FFT spatial derivatives with zero Nyquist create extra checkerboard null modes. Smooth manufactured modes avoid them; this discretization is not a production proposal.',
 'CG has no preconditioner. CPU cost is not a GPU or MPI estimate.',
 'No variable collisions, perpendicular magnetic driver, saturation, or pressure evolution.']
Path(__file__).with_name('rational_results.json').write_text(json.dumps(results,indent=2)+'\n')
assert results['cases'][4]['rms_error_relative_or_null_absolute'] < results['cases'][0]['rms_error_relative_or_null_absolute']
