/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    BioGrowth


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

#include "speciesTable.H"

#include "reactingWallFvPatchScalarField.C"
#include "reactingWallFvPatchScalarField.H"
using namespace std;

extern"C" {
void invokebrns_(double *theCurArray, double *thePreArray, double *outputArray, int *numComp, 
	double *time_step, int *boundary_flag, int *return_value, double *x_pos, double *y_pos, 
	double *z_pos, double *porosity, double *saturation, double *parameterVector);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #   include "setRootCase.H"

    #   include "createTime.H"
    #   include "createMesh.H"

    simpleControl simple(mesh);

    #   include "createFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

#   include "CourantNo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	double pos_x = 1;
	double pos_y = 1;
	double pos_z = 1;

	double porosity = 1;
	double waterSaturation = 1;

	int numComp = 4;

//	double concentrationTemp_;
	double* concentrationTemp_ = NULL;
	concentrationTemp_ = new double[numComp];

	int* boundary_flag  = NULL;
	boundary_flag = new int[numComp];

	int return_value   = -1;

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 	Info << "setting values of bio's Concentration around patch wall to 0.1." << endl;
	
	forAll(bio.boundaryField(), patchi)
	{
		if (Surf.boundaryField()[patchi].type()=="reactingWall")
		{
			Info << "boundary field type is: " << Surf.boundaryField()[patchi].type() << endl;
			
			const labelList& cellOwner = bio.boundaryField()[patchi].patch().faceCells();
			
			Info << "number of cells in cellowner: " << cellOwner.size() << endl;

			forAll(bio.boundaryField()[patchi], facei)
			{
				bio[cellOwner[facei]] = 0.1;
			}
		}
	}
	
	Info << "End of bio's Concentraion update." << endl;
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
		#include "CEqn.H"
        }


        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}



// ************************************************************************* //
