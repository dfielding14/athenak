from pathlib import Path
p=Path(__file__).parent/'source/src/pgen/tests/cgl_fofc.cpp'; s=p.read_text()
s=s.replace('  const bool weak_field = mode == "weak_field_transport";', '''  const bool affine = mode == "q_affine_compression";
  const bool weak_field = mode == "weak_field_transport" || affine;
  const int affine_dir = pin->GetOrAddInteger("problem", "q_affine_direction", 0);
  const Real affine_rate = pin->GetOrAddReal("problem", "q_affine_rate", -0.2);''')
s=s.replace('  auto b0 = pmhd->b0;\n\n  par_for("cgl_fofc_e2e_init"', '''  auto b0 = pmhd->b0;
  const Real affine_xmin = pmy_mesh_->mesh_size.x1min;
  const Real affine_dx = (pmy_mesh_->mesh_size.x1max-affine_xmin)/indcs.nx1;
  const Real affine_bx = affine_dir == 0 ? 10.0 : 0.0;
  const Real affine_by = affine_dir == 1 ? 10.0 : (affine_dir == 2 ? 1.0e-12 : 0.0);

  par_for("cgl_fofc_e2e_init"''')
s=s.replace('      bcc0(m,IBY,k,j,i) = weak_side ? 1.0e-12 : 1.0;\n    }', '''      bcc0(m,IBY,k,j,i) = weak_side ? 1.0e-12 : 1.0;
    }
    if (affine) {
      w0(m,IDN,k,j,i) = 1.0;
      w0(m,IVX,k,j,i) = affine_rate*(affine_xmin+(i-is+0.5)*affine_dx);
      w0(m,IVY,k,j,i) = w0(m,IVZ,k,j,i) = 0.0;
      w0(m,IPR,k,j,i) = w0(m,IPP,k,j,i) = 1.0;
      bcc0(m,IBX,k,j,i) = affine_bx;
      bcc0(m,IBY,k,j,i) = affine_by;
      bcc0(m,IBZ,k,j,i) = 0.0;
    }''')
s=s.replace('    b0.x1f(m,k,j,i) = weak_field ? 0.0 : (pressure_step ? 2.0 : 0.43);', '    b0.x1f(m,k,j,i) = affine ? affine_bx : (weak_field ? 0.0 : (pressure_step ? 2.0 : 0.43));')
s=s.replace('      b0.x2f(m,k,j,i) = weak_side ? 1.0e-12 : 1.0;', '      b0.x2f(m,k,j,i) = affine ? affine_by : (weak_side ? 1.0e-12 : 1.0);')
p.write_text(s)
