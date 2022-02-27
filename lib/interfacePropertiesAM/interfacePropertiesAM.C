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

#include "interfacePropertiesAM.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvCFD.H"
// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::interfacePropertiesAM::convertToRad =
    Foam::mathematicalConstant::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfacePropertiesAM::correctContactAngle
(
    surfaceVectorField::GeometricBoundaryField& nHatb,
    surfaceVectorField::GeometricBoundaryField& gradAlphaf
) const
{
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::GeometricBoundaryField& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleFvPatchScalarField& acap =
                const_cast<alphaContactAngleFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            vectorField nf = boundary[patchi].nf();

            // Reset nHatp to correspond to the contact angle

            scalarField a12 = nHatp & nf;

            scalarField b1 = cos(theta);

            scalarField b2(nHatp.size());

            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            scalarField det = 1.0 - a12*a12;

            scalarField a = (b1 - a12*b2)/det;
            scalarField b = (b2 - a12*b1)/det;

            nHatp = a*nf + b*nHatp;

            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}


void Foam::interfacePropertiesAM::calculateK()
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

	//init alpha smoothed
    volScalarField alpha1s = alpha1_;

	// smooth alpha1 by successive interpolation to face and cell
	// Raeini's thesis, equation (2.14)
	for (int i=0;i<nSK_;i++)
	{
    	alpha1s = cSK_ * fvc::average(fvc::interpolate(alpha1s)) + (1.0 - cSK_) * alpha1s;
	}
    
    // Cell gradient of alpha
	//volVectorField gradAlpha = fvc::grad(alpha1_);
    volVectorField gradAlpha = fvc::grad(alpha1s,"nHat");
    //update interface vector at cell center, used for filtering
	nI_ = gradAlpha/(Foam::mag(gradAlpha) + deltaN_);
	/*----------------------------------------------*/

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf = fvc::interpolate(gradAlpha);
    //gradAlphaf -=
    //    (mesh.Sf()/mesh.magSf())
    //   *(fvc::snGrad(alpha1_) - (mesh.Sf() & gradAlphaf)/mesh.magSf());

    // Face unit interface normal
    //surfaceVectorField nHatfv = gradAlphaf/(mag(gradAlphaf) + deltaN_);
    surfaceVectorField nHatfv = fvc::interpolate(nI_);
    correctContactAngle(nHatfv.boundaryField(), gradAlphaf.boundaryField());


    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    //K_ = -fvc::div(nHatf_);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    K_ = -fvc::div(nHatf_) + (nI_ & fvc::grad(nHatfv) & nI_);

    /*
    volVectorField nHat = gradAlpha/(mag(gradAlpha) + deltaN_);
    forAll(nHat.boundaryField(), patchi)
    {
        nHat.boundaryField()[patchi] = nHatfv.boundaryField()[patchi];
    }

    K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
    */
    
	// Curvature at the center of faces
	Kf_ = fvc::interpolate(K_);
}

//Re-calculate capillary flux
void Foam::interfacePropertiesAM::calculatePhic()
{
	const fvMesh& mesh = alpha1_.mesh();
	const surfaceScalarField& magSf = mesh.magSf();
 
    volScalarField alpha_pc
	(
		IOobject 
		(
			"alpha_pc",
			alpha1_.time().timeName(),
			alpha1_.mesh()
		),
	 alpha1_,
	 pc_.boundaryField().types()
	);

    // Sharpen interface function
	// Raeini's thesis, equation (2.22) and (2.23)
    alpha_pc = 1.0/(1.0-cPc_)*(min( max(alpha1_,cPc_/2.0), (1.0-cPc_/2.0) ) - cPc_/2.0);

	alpha_pc.correctBoundaryConditions();

    surfaceScalarField deltasf = fvc::snGrad(alpha_pc);

	//surface tension force
    surfaceScalarField stf = sigma_*Kf_*deltasf;

	//delta interface, used for filtering
	surfaceScalarField dInterf=(mag(deltasf)/(mag(deltasf)+1.0e-4*deltaN_));

	//surface tension force flux
	phic_ = stf*magSf -phicFilt1_*dInterf;

    for(int nonOrth=0; nonOrth<=nNonOrthogonalCorrectors_; nonOrth++)
    {
		//solve for pc
        fvScalarMatrix pcEqn
        (
            fvm::laplacian( pc_) == fvc::div(phic_)
        );

        pcEqn.setReference(pcRefCell_, pcRefValue_);

        pcEqn.solve();

		//add flux of pc to capillary flux
        if (nonOrth == nNonOrthogonalCorrectors_)
        { 
            phic_-=pcEqn.flux();
        }
    }

	//reconstruct capillary field   
    gPc_= fvc::reconstruct(phic_);
    gPc_.correctBoundaryConditions();


	//gradPc -(gradPc.nS)nS for equation (2.25), Raeini's thesis
	volVectorField vgPcllInterface=(gPc_)-((gPc_)&(nI_))*(nI_);

	//<gradPc -(gradPc.nS)nS>_f for equation (2.25), Raeini's thesis
	surfaceVectorField gPcllInterface=fvc::interpolate(vgPcllInterface);

	//Equation (2.25), Raeini's thesis
	phicFilt1_=dInterf*(cfcFiltRelax_*phicFilt1_+cfcFilt_*(gPcllInterface& mesh.Sf())*( 1.0-mag(nHatf_)/magSf));

	//if phic and phicFilt1 opposite sign, make sure |phicFilt1|<|phic|
	forAll(phicFilt1_,i)
	{
	    phicFilt1_[i]-=(phic_[i]*phicFilt1_[i] < 0) ? (max(min(phicFilt1_[i],mag(phic_[i])),-mag(phic_[i]))) : 0.0;
	}
	const fvBoundaryMesh& boundary = mesh.boundary(); 
	forAll(boundary, patchI)
	{
		scalarField phicFilt1bf = phicFilt1_.boundaryField()[patchI];
		scalarField phicbf = phic_.boundaryField()[patchI];
		forAll(phicFilt1bf,i)
		{
	        phicFilt1bf[i]-=(phicbf[i]*phicFilt1bf[i] < 0) ? (max(min(phicFilt1bf[i],mag(phicbf[i])),-mag(phicbf[i]))) : 0.0;
        }
	}

	//compute average of surface tension force on interface for second filtering
	dimensionedScalar stfAvg = dimensionedScalar("stfavg", dimPressure / dimLength, 0.0);
	stfAvg = average(mag(stf)*dInterf)/average(dInterf+1e-9);
	surfaceScalarField phicThreshold = cphicFilt_*(dInterf/(dInterf+1e-9))*stfAvg*magSf;
	//second filtering. Raeini's thesis, equation (2.26)
	phicFilt2_ = (min (max(phic_, -phicThreshold),phicThreshold));
	phic_-=phicFilt2_;

	//recompute capillary field
	gPc_= fvc::reconstruct(phic_);
	gPc_.correctBoundaryConditions();

	//additional capillary flux filtering
	if (gPcCorr_)
	{
		//smooth capillary field
		volVectorField gPcS=fvc::average(fvc::interpolate(gPc_,"limitedScheme"));
		
		//additional filtering, not in Raeini's thesis.
		phicThreshold = 1.5*mag((fvc::interpolate(gPcS)) & mesh.Sf());
		phic_ = (min (max(phic_, -phicThreshold),phicThreshold));

		gPc_= fvc::reconstruct(phic_);
		gPc_.correctBoundaryConditions();
	}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacePropertiesAM::interfacePropertiesAM
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cAlpha_
    (
        readScalar
        (
            alpha1.mesh().solutionDict().subDict("PISO").lookup("cAlpha")
        )
    ),
	//- smoothing coefficient for alphaS (generally 0.5)
	cSK_
    (
        readScalar
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("cSK")
        )
    ),
	//-number of smoothing cycle for alphaS
	nSK_
    (
        readLabel
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("nSK")
        )
    ),
	//constants for filtering capillary force parallel to interface
    cfcFilt_
    (
        readScalar
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("cfcFilt")
        )
    ),
    cfcFiltRelax_
    (
        readScalar
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("cfcFiltRelax")
        )
    ),
	//constant for filtering net capillary flux
	cphicFilt_
    (
        readScalar
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("cphicFilt")
        )
    ),
	//additional capillary flux filtering
	gPcCorr_
    (
        readBool
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("gPcCorr")
        )
    ),
	//- Sharp force coefficient (put 0.98-0.99 for static problems, 0.4-0.5 for dynamic)
	cPc_
    (
        readScalar
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("cPc")
        )
    ),
	//-number of non-orthogonal corrector loop
	nNonOrthogonalCorrectors_
    (
        readLabel
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("nNonOrthogonalCorrectors")
        )
    ),
    sigma_(dict.lookup("sigma")),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
    ),

    alpha1_(alpha1),
    U_(U),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    ),

    K_
    (
        IOobject
        (
            "K",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0)
    ),
	//intreface normal vector at cell center
	nI_
	(
		IOobject
		(
			"nI",
			alpha1_.time().timeName(),
			alpha1_.mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		alpha1_.mesh(),
		dimensionedVector("nI", dimless, vector(0.0, 0.0, 0.0))
	),
	//interface curvature at face center
    Kf_
    (
        IOobject
        (
            "Kf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("Kf", dimless/dimLength, 0.0)
    ),
	//capillary pressure
	pc_
	(
		IOobject
		(
			"pc",
			alpha1_.time().timeName(),
			alpha1_.mesh(),
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		alpha1_.mesh()
	),
	//Specific cell where the capillary pressure value is used as refference
	pcRefCell_
	(
		readLabel
		(
			alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("pcRefCell")
		)
	),
	//value of the capillary pressure in the specific cell
	pcRefValue_
	(
		readScalar
		(
			alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("pcRefValue")
		)
	),
	//capillary flux
	phic_
	(
		IOobject
		(
			"phic",
			alpha1_.time().timeName(),
			alpha1_.mesh()
		),
		alpha1_.mesh(),
		dimensionedScalar("phic", dimPressure / dimLength*dimArea, 0.0)
	),

	//capillary field
    gPc_
    (
        IOobject
        (
            "gPc",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector("gPc", dimPressure/dimLength, vector(0.0,0.0,0.0)),
        pc_.boundaryField().types()
    ),
	//First filtering function
    phicFilt1_
    (
        IOobject
        (
            "phicFilt1",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar("phicFilt1", dimPressure/dimLength*dimArea, 0.0)
    ),

	//second filtering function
    phicFilt2_
    (
        IOobject
        (
            "phicFilt2",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar("phicFilt2", dimPressure/dimLength*dimArea, 0.0)
    )
{
    setRefCell
    (
        pc_,
        alpha1_.mesh().solutionDict().subDict("PIMPLE"),
        pcRefCell_,
        pcRefValue_
    );

    calculateK();
}


// ************************************************************************* //
