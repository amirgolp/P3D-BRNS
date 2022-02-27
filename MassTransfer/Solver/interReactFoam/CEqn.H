

volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    DT1 * DT2 / (alpha1 * DT1 + (1 - alpha1) * DT2) 
);

{
    surfaceScalarField phiCiH = 
        (
            ( 
                fvc::interpolate(DT) * 
                (1-He)
               / (fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*He)
            )
            * fvc::snGrad(alpha1)
        ) * mesh.magSf();

    solve
    (
        fvm::ddt(H)
      + fvm::div(phi, H, "div(phi,C)")
      - fvm::laplacian(fvc::interpolate(DT), H, "laplacian(C)")
      + fvm::div(phiCiH, H, "div(phi,C)")
     ,
        mesh.solver("H")
    );
}

{
    surfaceScalarField phiCiM =
    (
        ( DT1 - DT2*He )
      / (fvc::interpolate(alpha1) + (1-fvc::interpolate(alpha1))*He)
      * fvc::snGrad(alpha1)
    ) * mesh.magSf();

    surfaceScalarField phiD = 
    (
        fvc::snGrad(DT) 
    ) * mesh.magSf();
    
    solve
    (
        fvm::ddt(M)
      + fvm::div(phi, M, "div(phi,C)")
      - fvm::laplacian(fvc::interpolate(DT), M, "laplacian(C)")
      - fvm::div(phiD, M, "div(phi,C)")
      + fvm::div(phiCiM, M, "div(phi,C)")
    ,
        mesh.solver("M")
    );
}
