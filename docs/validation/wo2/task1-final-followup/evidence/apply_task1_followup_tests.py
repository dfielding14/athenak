from pathlib import Path
import difflib,hashlib,json,ast
here=Path(__file__).resolve().parent
prod=Path('/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2')
changes={
 'tst/test_suite/cgl/test_cgl_lf_acceptance_cpu.py':[(
 '    result = run_case(tmp_path, name, "problem/validation_output=true",\n                      f"problem/validation_output_dir={tmp_path}")',
 '    # For this aligned, fixed-coefficient closure, dt_FE = dx**2/(2*chi_parallel).\n    # With sts_safety=0.9, a cap of 8 gives 2*128**2/(0.9*8*(2*pi)**2) > 115\n    # full-step equivalents per parallel damping time; retain >=100 resolved\n    # cycles without changing the physical duration or either closure coefficient.\n    result = run_case(tmp_path, name, "time/sts_max_dt_ratio=8",\n                      "problem/validation_output=true",\n                      f"problem/validation_output_dir={tmp_path}")')],
 'tst/test_suite/cgl/test_cgl_landau_fluid_cpu.py':[(
 '        history = testutils.athena_read.hst("cgl_ci_low_field.mhd.hst")\n        assert history["lf_nstage"][-1] > 0.0',
 '        history = testutils.athena_read.hst("cgl_ci_low_field.mhd.hst")\n        # Every face is disabled below bfloor, so the discrete LF row is zero.\n        # The timestep controller skips LF stages; the pgen still checks the\n        # physical pressure state against its initial values at 1e-13.\n        assert history["lf_nstage"][-1] == 0.0')],
 'tst/test_suite/cgl/test_cgl_landau_fluid_mpicpu.py':[(
 '        for heating_rate, post_stages in ((0, 7), (3000, 11)):',
 '        # Task1\'s aligned row gives dt_FE=h^2/(2*chi), chi=sqrt(8/pi)/(2*pi).\n        # The cap is 20*0.9*dt_FE; p_after=1+(2/3)*heating_rate*dt_cycle.\n        # A half-sweep ratio 10*sqrt(p_after) requires 7 or 13 RKL2 stages.\n        for heating_rate, post_stages in ((0, 7), (3000, 13)):')],
 'tst/test_suite/cgl/test_cgl_amr_gpu.py':[(
 '        _run("cgl_lf_amr_primitive_churn.athinput", "cgl_amr_gpu_lf_churn")',
 '        # Cap only the cycle duration so the fixed physical-time threshold\n        # still samples both refinement and derefinement with the new LF bound.\n        _run("cgl_lf_amr_primitive_churn.athinput", "cgl_amr_gpu_lf_churn",\n             "time/sts_max_dt_ratio=0.1")'),(
 '        _run("cgl_lf_amr_3d_current.athinput", "cgl_amr_gpu_lf_3d_churn")',
 '        # Preserve the original cycle-three AMR event at fixed switch time.\n        _run("cgl_lf_amr_3d_current.athinput", "cgl_amr_gpu_lf_3d_churn",\n             "time/sts_max_dt_ratio=0.13")')],
 'tst/test_suite/cgl/test_cgl_amr_walls_cpu.py':[(
 '         "time/nlim=3", "output2/variable=mhd_w_bcc", "output2/id=state"],',
 '         # Restore the original cycle-three event using only a timestep cap.\n         "time/sts_max_dt_ratio=0.13", "time/nlim=3",\n         "output2/variable=mhd_w_bcc", "output2/id=state"],')],
 'tst/test_suite/cgl/test_cgl_amr_mpi_gpu.py':[(
 '        _launch(2, 16, "cgl_lf_amr_3d_current.athinput", basename)',
 '        # Match the fixed physical AMR switch under Task1\'s larger LF bound.\n        _launch(2, 16, "cgl_lf_amr_3d_current.athinput", basename,\n                "time/sts_max_dt_ratio=0.13")')],
}
patch=[];manifest={'classification':'Task1 final regression adaptation; input decks and scientific assertions unchanged','files':{}}
for rel,replacements in changes.items():
 original=(here/'source-fused'/rel).read_text(); scratch=original
 production=(prod/rel).read_text()
 for old,new in replacements:
  assert scratch.count(old)==1,rel
  scratch=scratch.replace(old,new)
  if old in production:
   assert production.count(old)==1,rel
   production=production.replace(old,new)
  else:assert new in production,rel
 ast.parse(scratch,filename=rel)
 (prod/rel).write_text(production)
 (here/'regression-source'/rel).write_text(scratch)
 patch.extend(difflib.unified_diff(original.splitlines(True),scratch.splitlines(True),fromfile='a/'+rel,tofile='b/'+rel))
 manifest['files'][rel]={'before_sha256':hashlib.sha256(original.encode()).hexdigest(),'after_sha256':hashlib.sha256(scratch.encode()).hexdigest()}
(here/'task1-final-regression.patch').write_text(''.join(patch))
manifest['patch_sha256']=hashlib.sha256((here/'task1-final-regression.patch').read_bytes()).hexdigest()
(here/'task1-final-regression.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(manifest['patch_sha256'])
