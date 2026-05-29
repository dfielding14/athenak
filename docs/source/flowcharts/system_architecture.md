# AthenaK System Architecture

The diagrams below mirror the control flow in `src/main.cpp`, `src/mesh/meshblock_pack.cpp`, and the physics constructors. Nodes link directly to the relevant documentation pages so you can drill into the implementation.

## Layered Overview

```{mermaid}
flowchart LR
    classDef input fill:#e3f2fd,stroke:#1e88e5,stroke-width:2px,color:#0d47a1;
    classDef boot fill:#fff3e0,stroke:#ef6c00,stroke-width:2px,color:#e65100;
    classDef pack fill:#e8f5e9,stroke:#388e3c,stroke-width:2px,color:#1b5e20;
    classDef physics fill:#f3e5f5,stroke:#6a1b9a,stroke-width:2px,color:#4a148c;
    classDef support fill:#fce4ec,stroke:#ad1457,stroke-width:2px,color:#880e4f;
    classDef output fill:#ede7f6,stroke:#4527a0,stroke-width:2px,color:#311b92;

    subgraph INPUT["Input sources"]
        CLI["CLI flags (-i, -r, -t, ...)"]
        ATHINPUT["athinput file"]
        RESTART["Restart dataset"]
        XDATA["Problem-specific data"]
    end
    class CLI,ATHINPUT,RESTART,XDATA input

    subgraph BOOT["Boot sequence"]
        MAIN["main.cpp"]
        PARAMS["ParameterInput"]
        MESH["Mesh & MeshBlockTree"]
    end
    class MAIN,PARAMS,MESH boot

    subgraph PACKLAYER["Mesh::AddCoordinatesAndPhysics"]
        COORDS["Coordinates"]
        PACK["MeshBlockPack"]
    end
    class PACK,COORDS pack

    subgraph PHYS["Physics modules"]
        HYDRO["Hydro"]
        MHD["MHD"]
        ION["Ion-neutral\n(requires Hydro + MHD)"]
        RAD["Radiation"]
        Z4C["Z4c"]
        ADM["ADM"]
        GRMHD["Dyn GRMHD"]
        PARTICLES["Particles"]
        TURB["Turbulence driver"]
    end
    class HYDRO,MHD,ION,RAD,Z4C,ADM,GRMHD,PARTICLES,TURB physics

    subgraph SUPPORTS["Shared services"]
        SRCTERMS["Source terms"]
        DIFF["Diffusion\n(viscosity / resistivity / conduction)"]
        BVALS["Boundary exchange"]
        TASKS["TaskList map"]
        PGEN["Problem generator"]
    end
    class SRCTERMS,DIFF,BVALS,TASKS,PGEN support

    subgraph OUTPUTLAYER["Outputs"]
        DRIVER["Driver"]
        OUTPUTS["Outputs manager"]
        HIST["History"]
        RST["Restart writer"]
        VTK["VTK / pvtk"]
        BIN["Binary / coarsened"]
    end
    class DRIVER,OUTPUTS,HIST,RST,VTK,BIN output

    CLI --> MAIN
    ATHINPUT --> PARAMS
    RESTART --> PARAMS
    XDATA --> PGEN

    MAIN --> PARAMS
    PARAMS --> MESH
    MESH --> PACK

    PACK --> COORDS
    PACK --> HYDRO
    PACK --> MHD
    PACK --> ION
    PACK --> RAD
    PACK --> Z4C
    PACK --> ADM
    PACK --> GRMHD
    PACK --> PARTICLES
    PACK --> TURB

    HYDRO --> SRCTERMS
    HYDRO --> DIFF
    MHD --> SRCTERMS
    MHD --> DIFF
    RAD --> TASKS
    TURB --> TASKS
    Z4C --> GRMHD
    ADM --> GRMHD

    MESH --> PGEN
    PGEN --> HYDRO
    PGEN --> MHD
    PGEN --> RAD
    PGEN --> PARTICLES

    HYDRO --> TASKS
    MHD --> TASKS
    ION --> TASKS
    RAD --> TASKS
    GRMHD --> TASKS
    PARTICLES --> TASKS
    Z4C --> TASKS
    ADM --> TASKS

    TASKS --> DRIVER
    DRIVER --> OUTPUTS
    OUTPUTS --> HIST
    OUTPUTS --> RST
    OUTPUTS --> VTK
    OUTPUTS --> BIN

    click HYDRO "../modules/hydro.html" "Hydrodynamics module"
    click MHD "../modules/mhd.html" "MHD module"
    click ION "../modules/ion_neutral.html" "Ion-neutral driver"
    click RAD "../modules/radiation.html" "Radiation module"
    click Z4C "../modules/z4c.html" "Z4c metric"
    click GRMHD "../modules/dyn_grmhd.html" "Dynamical GRMHD"
    click PARTICLES "../modules/particles.html" "Particles module"
    click TURB "../modules/srcterms.html#turbulence-driver" "Turbulence driver"
    click COORDS "../modules/coordinates.html" "Coordinates system"
    click TASKS "../modules/tasklist.html" "TaskList system"
    click SRCTERMS "../modules/srcterms.html" "Source terms overview"
    click DIFF "../modules/diffusion.html" "Diffusion and resistivity"
    click BVALS "../modules/boundaries.html" "Boundary exchange"
    click PGEN "../modules/pgen.html" "Problem generators"
    click DRIVER "../modules/driver.html" "Driver lifecycle"
    click OUTPUTS "../modules/outputs.html" "Outputs manager"
```

Ion-neutral MHD is instantiated only when `<ion-neutral>`, `<hydro>`, and
`<mhd>` blocks all exist. Specifying both fluid blocks without
`<ion-neutral>` also aborts. Dynamical spacetime inputs (`<z4c>` or `<adm>`)
are rejected with Hydro and with ion-neutral coupling, while an MHD state
routes through DynGRMHD.

## Task Contributors

The `TaskList` map inside each `MeshBlockPack` aggregates stage tasks from every active module. The diagram below highlights the main contributors and their dependencies.

```{mermaid}
flowchart LR
    classDef core fill:#fff3e0,stroke:#ef6c00,stroke-width:2px,color:#bf360c;
    classDef phys fill:#f3e5f5,stroke:#6a1b9a,stroke-width:2px,color:#4a148c;
    classDef support fill:#fce4ec,stroke:#ad1457,stroke-width:2px,color:#880e4f;

    TL["TaskList scheduler\n(before/stage/after)"]:::core
    HYDRO_TASKS["Hydro tasks\n(AssembleHydroTasks)"]:::phys
    MHD_TASKS["MHD tasks\n(AssembleMHDTasks)"]:::phys
    ION_TASKS["Ion-neutral tasks\n(AssembleIonNeutralTasks)"]:::phys
    RAD_TASKS["Radiation tasks"]:::phys
    NUMREL["Numerical relativity tasks"]:::phys
    PART["Particle tasks"]:::phys
    TURB_TASKS["Turbulence driver tasks\n(InitializeModes, AddForcing)"]:::phys

    BVAL_TASKS["Boundary exchange tasks"]:::support
    DIFF_TASKS["Diffusion / resistivity"]:::support
    SRC_TASKS["Source terms (Hydro/MHD)"]:::support

    TL --> HYDRO_TASKS
    TL --> MHD_TASKS
    TL --> ION_TASKS
    TL --> RAD_TASKS
    TL --> NUMREL
    TL --> PART
    TL --> TURB_TASKS
    HYDRO_TASKS --> SRC_TASKS
    HYDRO_TASKS --> BVAL_TASKS
    HYDRO_TASKS --> DIFF_TASKS
    MHD_TASKS --> SRC_TASKS
    MHD_TASKS --> BVAL_TASKS
    MHD_TASKS --> DIFF_TASKS
    ION_TASKS --> MHD_TASKS
    NUMREL --> MHD_TASKS
    RAD_TASKS --> BVAL_TASKS
    PART --> BVAL_TASKS
    TURB_TASKS --> HYDRO_TASKS
    TURB_TASKS --> MHD_TASKS

    class HYDRO_TASKS,MHD_TASKS,ION_TASKS,RAD_TASKS,NUMREL,PART,TURB_TASKS phys
    class BVAL_TASKS,DIFF_TASKS,SRC_TASKS support

    click TL "../modules/tasklist.html"
    click HYDRO_TASKS "../modules/hydro.html"
    click MHD_TASKS "../modules/mhd.html"
    click ION_TASKS "../modules/ion_neutral.html"
    click RAD_TASKS "../modules/radiation.html"
    click NUMREL "../modules/dyn_grmhd.html"
    click PART "../modules/particles.html"
    click TURB_TASKS "../modules/srcterms.html#turbulence-driver"
    click BVAL_TASKS "../modules/boundaries.html"
    click DIFF_TASKS "../modules/diffusion.html"
    click SRC_TASKS "../modules/srcterms.html"
```

## MeshBlock Data Flow

The following diagram summarises the major arrays owned by each `MeshBlockPack` and how they feed the physics kernels. It reflects the data layout used by Hydro/MHD and the boundary system (`src/hydro/hydro.cpp`, `src/mhd/mhd.cpp`, `src/bvals`).

```{mermaid}
flowchart LR
    classDef data fill:#fffde7,stroke:#f9a825,stroke-width:2px,color:#f57f17;
    classDef compute fill:#c8e6c9,stroke:#388e3c,stroke-width:2px,color:#1b5e20;
    classDef comm fill:#d1c4e9,stroke:#4527a0,stroke-width:2px,color:#311b92;

    subgraph PACKDATA[Per MeshBlockPack]
        U0["u0 (conserved)"]
        W0["w0 (primitive)"]
        BFACE["b0.{x1f,x2f,x3f}"]
        BCC["bcc0"]
        SCALARS["Passive scalars"]
    end
    class U0,W0,BFACE,BCC,SCALARS data

    subgraph COMPUTE[Local computation]
        C2P[Conservative → Primitive]
        RECON[Reconstruction]
        RIEMANN[Riemann solvers]
        UPDATE[RK update]
        CT[Constrained transport]
    end
    class C2P,RECON,RIEMANN,UPDATE,CT compute

    subgraph COMM[MPI / boundary exchange]
        PACKBUF["Pack buffers"]
        EXCHANGE["Send and receive"]
        UNPACK["Unpack buffers"]
    end
    class PACKBUF,EXCHANGE,UNPACK comm

    U0 --> C2P
    C2P --> W0
    W0 --> RECON
    RECON --> RIEMANN
    RIEMANN --> UPDATE
    UPDATE --> U0
    BFACE --> CT
    CT --> BFACE
    CT --> BCC
    BCC --> RECON

    U0 --> PACKBUF
    BFACE --> PACKBUF
    PACKBUF --> EXCHANGE
    EXCHANGE --> UNPACK
    UNPACK --> U0
    UNPACK --> BFACE

    click RECON "../modules/reconstruction.html"
    click RIEMANN "../modules/riemann_solvers.html"
    click CT "../modules/mhd.html#constrained-transport"
    click EXCHANGE "../modules/boundaries.html"
```

These diagrams collectively describe how the bootstrap code wires physics modules into the task scheduler and how state flows during integration. For the temporal evolution of those tasks, see the runtime flowchart.
