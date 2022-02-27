/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef interfacePropertiesAM_H
#define interfacePropertiesAM_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interfacePropertiesAM Declaration
\*---------------------------------------------------------------------------*/

class interfacePropertiesAM
{

private:
    // Private data

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;

        //- Compression coefficient
        scalar cAlpha_;

		//- smoothing coefficient for alphaS
		scalar cSK_;
		//-number of smoothing cycle for alphaS
		label nSK_;
		//-Filtering constants
		scalar cfcFilt_;
		scalar cfcFiltRelax_;
        scalar cphicFilt_;
		bool   gPcCorr_;
		//- Sharp force coefficient (put 0.98-0.99 for static problems, 0.4-0.5 for dynamic)
		scalar cPc_;
		//-number of non-orthogonal corrector loop
		label nNonOrthogonalCorrectors_;
		/*----------------------------------------------*/

        //- Surface tension
        dimensionedScalar sigma_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        const volScalarField& alpha1_;
        const volVectorField& U_;
        surfaceScalarField nHatf_;
        volScalarField K_;

		//interface normal vector at cell center
		volVectorField nI_;
		//interface curvature at face center
        surfaceScalarField Kf_;
		//capillary pressure
		volScalarField pc_;
		label pcRefCell_;
		scalar pcRefValue_;
		//Capillary flow	
		surfaceScalarField phic_;
		//capillary field gPc_ = stf-gradPc
		volVectorField gPc_;
		//functions for filtering capillary force parallel to interface
		surfaceScalarField phicFilt1_;
		surfaceScalarField phicFilt2_;
		/*----------------------------------------------*/

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        interfacePropertiesAM(const interfacePropertiesAM&);
        void operator=(const interfacePropertiesAM&);

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
            surfaceVectorField::GeometricBoundaryField& nHat,
            surfaceVectorField::GeometricBoundaryField& gradAlphaf
        ) const;

        //- Re-calculate the interface curvature
        void calculateK();

		//Re-calculate capillary flux
		void calculatePhic();

public:

    //- Conversion factor for degrees into radians
    static const scalar convertToRad;


    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        interfacePropertiesAM
        (
            const volScalarField& alpha1,
            const volVectorField& U,
            const IOdictionary&
        );


    // Member Functions

        scalar cAlpha() const
        {
            return cAlpha_;
        }

        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }

        const surfaceScalarField& nHatf() const
        {
            return nHatf_;
        }

        const volScalarField& K() const
        {
            return K_;
        }

        const dimensionedScalar& sigma() const
        {
            return sigma_;
        }

        tmp<volScalarField> sigmaK() const
        {
            return sigma_*K_;
        }

		//function correct now include calculate capillary flux
        void correct()
        {
            calculateK();
            calculatePhic();
        }

		//get capillary pressure
		const volScalarField& pc()
		{
			return pc_;
		}

		//get gradient capillary pressure
		const volVectorField& gPc()
		{
			return gPc_;
		}
		//get capillary flow of pc
		const surfaceScalarField& phic()
		{
			return phic_;
		}
		/*----------------------------------------------*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
