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

 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//if (std::abs(curTime - reactionDeltaT) < 1e-5 )
//{
	//reactionDeltaT += 1e-4;
	Info << "Reaction Step" << endl;

	double celli = 0;
	
	forAll(Y[0].boundaryField(), patchi)
	{
		
		if (Surf.boundaryField()[patchi].type()=="reactingWall")
		{

			const labelList& cellOwner = Y[0].boundaryField()[patchi].patch().faceCells();
			
			forAll(Y[0].boundaryField()[patchi], facei)
			{
				
				celli = cellOwner[facei];
				double time_step = runTime.deltaT().value();

				/*concentrationTemp_[0] = Y[0][celli];
				concentrationTemp_[1] = Y[1][celli];
				concentrationTemp_[2] = Y[2][celli];				
				concentrationTemp_[3] = Y[3][celli];				
				concentrationTemp_[4] = Y[4][celli];
				concentrationTemp_[5] = Y[5][celli];
				concentrationTemp_[6] = bio[celli];
				concentrationTemp_[7] = EsatNAP[celli];
				concentrationTemp_[8] = EsatNIR[celli];
				concentrationTemp_[9] = EsatNOR[celli];*/
				
				concentrationTemp_[0] = (Y[0][celli]>1e-20 ? Y[0][celli] : 0);
				concentrationTemp_[1] = (Y[1][celli]>1e-20 ? Y[1][celli] : 0);
				concentrationTemp_[2] = (Y[2][celli]>1e-20 ? Y[2][celli] : 0);
				concentrationTemp_[3] = (Y[3][celli]>1e-20 ? Y[3][celli] : 0);
				concentrationTemp_[4] = (Y[4][celli]>1e-20 ? Y[4][celli] : 0);
				concentrationTemp_[5] = (Y[5][celli]>1e-20 ? Y[5][celli] : 0);
				concentrationTemp_[6] = (bio[celli]>1e-20 ? bio[celli] : 0);
				concentrationTemp_[7] = (EsatNAP[celli]>1e-20 ? EsatNAP[celli] : 0);
				concentrationTemp_[8] = (EsatNIR[celli]>1e-20 ? EsatNIR[celli] : 0);
				concentrationTemp_[9] = (EsatNOR[celli]>1e-20 ? EsatNOR[celli] : 0);

				for (int i = 0 ; i < numComp; i++)
				{
					boundary_flag[i] = 0;
				}

		// Info << "Invoking BRNS" << endl;
				if (alpha1[celli] > 0.5) {
				invokebrns_
				(
					concentrationTemp_,
					concentrationTemp_,
					concentrationTemp_,
					&numComp,
					&time_step,
					boundary_flag,
					&return_value,
					&pos_x,
					&pos_y,
					&pos_z,
					&porosity,
					&waterSaturation,
					NULL
				);
				}
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
				/*Y[0][celli]    = (concentrationTemp_[0]);//mol/m3
				Y[1][celli]    = (concentrationTemp_[1]);//mol/m3
				Y[2][celli]    = (concentrationTemp_[2]);//mol/m3
				Y[3][celli]    = (concentrationTemp_[3]);//mol/m3
				Y[4][celli]    = (concentrationTemp_[4]);//mol/m3
				Y[5][celli]    = (concentrationTemp_[5]);//mol/m3
				bio[celli]     = (concentrationTemp_[6]);//mol/m3
				EsatNAP[celli]    = (concentrationTemp_[7]);//mol/m3
				EsatNIR[celli]    = (concentrationTemp_[8]);//mol/m3
				EsatNOR[celli]    = (concentrationTemp_[9]);//mol/m3*/
 
				Y[0][celli]    = (concentrationTemp_[0]>1e-20 ? concentrationTemp_[0] : 0);//mol/m3
				Y[1][celli]    = (concentrationTemp_[1]>1e-20 ? concentrationTemp_[1] : 0);//mol/m3
				Y[2][celli]    = (concentrationTemp_[2]>1e-20 ? concentrationTemp_[2] : 0);//mol/m3
				Y[3][celli]    = (concentrationTemp_[3]>1e-20 ? concentrationTemp_[3] : 0);//mol/m3
				Y[4][celli]    = (concentrationTemp_[4]>1e-20 ? concentrationTemp_[4] : 0);//mol/m3
				Y[5][celli]    = (concentrationTemp_[5]>1e-20 ? concentrationTemp_[5] : 0);//mol/m3
				bio[celli]     = (concentrationTemp_[6]>1e-20 ? concentrationTemp_[6] : 0);//mol/m3
				EsatNAP[celli]    = (concentrationTemp_[7]>1e-20 ? concentrationTemp_[7] : 0);//mol/m3
				EsatNIR[celli]    = (concentrationTemp_[8]>1e-20 ? concentrationTemp_[8] : 0);//mol/m3
				EsatNOR[celli]    = (concentrationTemp_[9]>1e-20 ? concentrationTemp_[9] : 0);//mol/m3
				
			}
		}
	}
//}	
}

Info << "End of reaction step" << endl;
Info << " " << endl;
