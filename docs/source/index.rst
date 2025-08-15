.. AthenaK documentation master file

============================================
AthenaK: Astrophysical Simulation Framework
============================================

.. raw:: html

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 2em; border-radius: 10px; margin: 1em 0;">
   <h2 style="color: white; margin: 0;">Performance-Portable Astrophysics on CPUs & GPUs</h2>
   <p style="color: #f0f0f0; margin-top: 0.5em;">Hydrodynamics • MHD • General Relativity • Radiation • Particles</p>
   </div>

System Overview - Start Here!
=============================

**Understand the complete system flow before diving into details:**

.. toctree::
   :maxdepth: 1
   
   overview

Quick Navigation
================

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 1em; margin: 2em 0;">
   
   <div style="border: 2px solid #667eea; padding: 1em; border-radius: 8px;">
   <h3>🚀 Getting Started</h3>
   <ul>
   <li><a href="quickstart.html">5-Minute Quickstart</a></li>
   <li><a href="building.html">Building AthenaK</a></li>
   <li><a href="configuration.html">Configuration Guide</a></li>
   <li><a href="running.html">Running Simulations</a></li>
   </ul>
   </div>
   
   <div style="border: 2px solid #764ba2; padding: 1em; border-radius: 8px;">
   <h3>📚 Complete Module Reference</h3>
   <ul>
   <li><a href="modules/index.html">All 20 Modules</a></li>
   <li><a href="reference/input_parameters.html">340 Input Parameters</a></li>
   <li><a href="reference/file_reference.html">File Reference</a></li>
   </ul>
   </div>
   
   <div style="border: 2px solid #667eea; padding: 1em; border-radius: 8px;">
   <h3>🔬 Physics Capabilities</h3>
   <ul>
   <li><a href="modules/hydro.html">Hydrodynamics</a></li>
   <li><a href="modules/mhd.html">Magnetohydrodynamics</a></li>
   <li><a href="modules/radiation.html">Radiation Transport</a></li>
   <li><a href="modules/z4c.html">Numerical Relativity</a></li>
   </ul>
   </div>
   
   <div style="border: 2px solid #764ba2; padding: 1em; border-radius: 8px;">
   <h3>⚡ Performance</h3>
   <ul>
   <li><a href="modules/tasklist.html">Task-Based Execution</a></li>
   <li><a href="modules/mesh.html">Adaptive Mesh Refinement</a></li>
   <li><a href="modules/boundaries.html">MPI Parallelization</a></li>
   </ul>
   </div>
   
   </div>

Complete Documentation
======================

.. toctree::
   :maxdepth: 2
   :caption: 🏁 Getting Started
   
   quickstart
   building
   configuration
   running
   troubleshooting

.. toctree::
   :maxdepth: 2
   :caption: 🔧 System Architecture
   
   flowcharts/system_architecture
   flowcharts/runtime
   kokkos_guide
   
.. toctree::
   :maxdepth: 1
   :caption: 📦 All Modules (20)
   
   modules/index
   modules/mesh
   modules/driver
   modules/tasklist
   modules/coordinates
   modules/hydro
   modules/mhd
   modules/radiation
   modules/z4c
   modules/dyn_grmhd
   modules/ion_neutral
   modules/particles
   modules/reconstruction
   modules/riemann_solvers
   modules/eos
   modules/diffusion
   modules/outputs
   modules/boundaries
   modules/srcterms
   modules/shearing_box
   modules/pgen

.. toctree::
   :maxdepth: 2
   :caption: 📖 Reference
   
   reference/input_parameters
   reference/file_reference
   reference/api_reference
   glossary
   
.. toctree::
   :maxdepth: 2
   :caption: 🔄 Migration
   
   migration/index
   migration/common_gotchas

.. toctree::
   :maxdepth: 1
   :caption: 📚 Examples
   
   examples/shock_tube
   examples/blast_wave
   examples/turbulence
   examples/mri_turbulence
   examples/binary_merger
   cgm_cooling_flow_metals

.. toctree::
   :maxdepth: 1
   :caption: 📝 Contributing
   
   contributing_docs

Key Features
============

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1em; margin: 2em 0;">
   
   <div style="text-align: center; padding: 1em;">
   <div style="font-size: 3em;">🚀</div>
   <h4>Performance Portable</h4>
   <p>Single code for CPUs & GPUs via Kokkos</p>
   </div>
   
   <div style="text-align: center; padding: 1em;">
   <div style="font-size: 3em;">🎯</div>
   <h4>Adaptive Mesh</h4>
   <p>Dynamic refinement for multi-scale physics</p>
   </div>
   
   <div style="text-align: center; padding: 1em;">
   <div style="font-size: 3em;">⚛️</div>
   <h4>Complete Physics</h4>
   <p>Hydro, MHD, GR, radiation, particles</p>
   </div>
   
   <div style="text-align: center; padding: 1em;">
   <div style="font-size: 3em;">⚡</div>
   <h4>Task-Based</h4>
   <p>Automatic dependency management</p>
   </div>
   
   </div>

Search & Index
==============

* :ref:`genindex`
* :ref:`search`

.. raw:: html

   <div style="margin-top: 3em; padding: 1em; background: #f5f5f5; border-left: 4px solid #667eea;">
   <strong>Need Help?</strong> Start with the <a href="overview.html">System Overview</a> to understand the architecture, 
   then explore specific <a href="modules/index.html">modules</a> as needed.
   </div>