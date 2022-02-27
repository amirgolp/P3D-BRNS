{

	//Transport step
	Info << "transport step" << endl;
	word alpharScheme("div(phirb,alpha)");
	word YiScheme("div(phi,Yi)");

    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(cYi*phic, max(phic));
    surfaceScalarField phir = phic*interface.nHatf();

    surfaceScalarField fluxDir = fvc::snGrad(alpha1)*mesh.magSf();
	surfaceScalarField alphaUp = upwind<scalar>(mesh,fluxDir).interpolate(alpha1);
	surfaceScalarField alphaDown = downwind<scalar>(mesh,fluxDir).interpolate(alpha1);

	//Re-initialize transfer flux
	phiD *=0;

    forAll(solutionSpecies, i)
	{
		volScalarField& Yi = Y[i];
		dimensionedScalar D1Yi = D1Y[i];
		dimensionedScalar D2Yi = D2Y[i];
		dimensionedScalar HYi = HY[i];

		Yi.correctBoundaryConditions();

		//compute diffusion coefficient
		surfaceScalarField D = fvc::interpolate(D1Yi*alpha1+HYi*D2Yi*(1-alpha1))/fvc::interpolate(alpha1+HYi*(1-alpha1));

		surfaceScalarField phiH = -D*(1-HYi)/fvc::interpolate((alpha1+HYi*(1-alpha1)))*fvc::snGrad(alpha1)*mesh.magSf();

		surfaceScalarField phiComp = fvc::flux
        (
            -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
            alpha1,
            alpharScheme
        );

		surfaceScalarField YiUp   = upwind<scalar>(mesh,fluxDir).interpolate(Yi);
		surfaceScalarField YiDown = downwind<scalar>(mesh,fluxDir).interpolate(Yi);
        surfaceScalarField deltaYi1 = max(-max(Yi),min(max(Yi),(YiDown-YiUp)/(alphaDown-alphaUp+1e-4)));
		surfaceScalarField deltaYi2 = YiDown/(alphaDown+(1-alphaDown)*HYi)-HYi*YiUp/(alphaUp+(1-alphaUp)*HYi);

		//interface compression flux
		surfaceScalarField phiY =
        fvc::flux
        (
            phi,
            Yi,
            YiScheme
        )
		+phiComp*sign(deltaYi2)*max(mag(deltaYi1),mag(deltaYi2));

    	MULES::explicitSolve(Yi, phi, phiY, maxY[i], 0);

		volScalarField SpY = (Yi-Yi.oldTime())/runTime.deltaT();

		while (simple.correctNonOrthogonal())
		{
			//compute Y equation
			fvScalarMatrix YiEqn
			(
				fvm::ddt(Yi)
		        - fvm::laplacian(D,Yi)
				- fvm::div(phiH,Yi,"div(phiH,Yi)")
				==
				SpY
			);

			//solve equations
		    if (simple.finalNonOrthogonalIter())
		    {
		        YiEqn.solve(mesh.solutionDict().solver("YiFinal"));
		    }
		    else
		    {
				YiEqn.solve(mesh.solutionDict().solver("Yi"));
			}
		}

		Info<< "Species concentration = "
        << Yi.weightedAverage(mesh.V()).value()
			<< "  Min(Yi) = " << gMin(Yi.internalField())
			<< "  Max(Yi) = " << gMax(Yi.internalField())
			<< endl;

		phiD+=Mw[i]*
			 (
			 	D*fvc::snGrad(Yi)*mesh.magSf()
			   +fvc::flux(phiH,Yi,"div(phiH,Yi)")
			 );
	}

	Mflux = fvc::div(phiD*alphaDown)-alpha1*fvc::div(phiD);

}

Info << "End of reaction step" << endl;
Info << " " << endl;
