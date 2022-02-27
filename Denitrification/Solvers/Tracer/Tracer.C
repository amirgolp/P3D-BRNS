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
    MoDesFinal

Description
    Steady-state solver for incompressible, turbulent flow, using the SIMPLE
    algorithm.

\*---------------------------------------------------------------------------*/
#include <cmath>
#include "fvCFD.H"
#include "simpleControl.H"
#include "pimpleControl.H"
#include "speciesTable.H"
#include "upwind.H"
#include "downwind.H"
#include "MULES.H"
#include "twoPhaseMixture.H"
//#include "write.H"

#include "reactingWallFvPatchScalarField.C"
#include "reactingWallFvPatchScalarField.H"
// #include "wallGlobalCompositionFvPatchScalarField.C"
#include "wallGlobalCompositionFvPatchScalarField.H"
//#include "interfacePropertiesAR.C"
//#include "interfacePropertiesAR.H"
#include "interfaceProperties.H"

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
    pimpleControl pimple(mesh);
    #   include "createFields.H"
    // #   include "readGravitationalAcceleration.H"
    #   include "createTimeControls.H"
    #   include "CourantNo.H"
    #   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



//double curTime = runTime.value();
//double reactionDeltaT = curTime + 1e-4;
//Info<< "Reaction deltaT is = " << reactionDeltaT << nl << endl;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"
	//curTime = runTime.value();

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (pimple.loop())
        {
	    twoPhaseProperties.correct();
			
            #include "CEqn.H"
        }
		
	runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}



// ************************************************************************* //
