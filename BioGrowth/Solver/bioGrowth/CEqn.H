

{
    fvScalarMatrix AEqn
    (
	fvm::ddt(A)
	+ fvm::div(phi, A, "div(phi,C)")
	- fvm::laplacian(DT, A, "laplacian(C)")
    );

    AEqn.solve();

    fvScalarMatrix BEqn
    (
	fvm::ddt(B)
	+ fvm::div(phi, B, "div(phi,C)")
	- fvm::laplacian(DT, B, "laplacian(C)")
    );

    BEqn.solve();

	fvScalarMatrix CEqn
    (
	fvm::ddt(C)
	+ fvm::div(phi, C, "div(phi,C)")
	- fvm::laplacian(DT, C, "laplacian(C)")
    );

    CEqn.solve();

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	double celli = 0;
	
	forAll(A.boundaryField(), patchi)
	{
		if (Surf.boundaryField()[patchi].type()=="reactingWall")
		{	
			const labelList& cellOwner = A.boundaryField()[patchi].patch().faceCells();
			
			forAll(A.boundaryField()[patchi], facei)
			{
				
				double time_step = runTime.deltaT().value();
				celli = cellOwner[facei];
				
				//set concentration for BRNS module
				//for (int i = 0 ; i < numComp; i++){
				concentrationTemp_[0] = A[celli];
				concentrationTemp_[1] = B[celli];
				concentrationTemp_[2] = 0;
				concentrationTemp_[3] = bio[celli];

				boundary_flag[0] = 0;
				boundary_flag[1] = 0;
				boundary_flag[2] = 0;
				boundary_flag[3] = 0;

				/*cout << "concentration A = " << concentrationTemp_[0] << endl;
				cout << "concentration B = " << concentrationTemp_[1] << endl;
				cout << "concentration C = " << concentrationTemp_[2] << endl;
				cout << "concentration bio = " << concentrationTemp_[3] << endl;*/

				// Run BRNS;
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

				if (concentrationTemp_[2] > 1e-20)
				{ 
					C[celli] += concentrationTemp_[2];
				}
				else if ( C[celli] > 1e-20)
				{
					C[celli] += 0;
				}
				else
				{
					C[celli] = 0;
				}

				A[celli]    = concentrationTemp_[0];//mol/m3
				B[celli]    = concentrationTemp_[1];//mol/m3
				bio[celli]  = concentrationTemp_[3];//mol/m3

			}
		}
	}		
}

