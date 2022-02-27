
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field alpha1\n" << endl;
volScalarField alpha1
(
    IOobject
    (
        "alpha1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "Reading/calculating face flux field phi\n" << endl;
#include "createPhi.H"
/*surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);*/
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

IOdictionary transportProperties
(
    IOobject
    (
	    "transportProperties",
	    runTime.constant(),
	    mesh,
	    IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*Info << "Reading component Succinate" << endl;
volScalarField C_s
(
    IOobject
    (
	"C_s",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading component Oxygen" << endl;
volScalarField C_O2
(
    IOobject
    (
	"C_O2",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
	//alpha1*scalar(0) + (scalar(1) - alpha1)*0.1
);

//forAll(C_O2, i) {
//		C_O2[i] = alpha1[i]*scalar(0) + (scalar(1) - alpha1[i])*0.1;	
//	}
//C_O2.write();

Info << "Reading component Nitrate" << endl;
volScalarField C_NO3
(
    IOobject
    (
	"C_NO3",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading component Nitrite" << endl;
volScalarField C_NO2
(
    IOobject
    (
	"C_NO2",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading component Nitric Oxide" << endl;
volScalarField C_NO
(
    IOobject
    (
	"C_NO",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading component Nitrous Oxide" << endl;
volScalarField C_N2O
(
    IOobject
    (
	"C_N2O",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);*/

Info << "Reading component Agrobacterium tumefaciens" << endl;
volScalarField bio
(
    IOobject
    (
	"bio",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading component Enzymes" << endl;
volScalarField EsatNAP
(
    IOobject
    (
	"EsatNAP",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField EsatNIR
(
    IOobject
    (
	"EsatNIR",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField EsatNOR
(
    IOobject
    (
	"EsatNOR",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

surfaceScalarField phiD
	(
		IOobject
		(
			"phiD",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("phiD", dimMass/dimTime, 0.0)
	); 

	volScalarField Mflux
	(
		IOobject
		(
			"Mflux",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("MFlux", dimMass/dimVolume/dimTime, 0.0)
	);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info << "Reading surface masters dictionary" << endl;
wordList surfaceMastersNames;
if (transportProperties.found("surfaceMasters"))
{
    const dictionary& surfacesDict = transportProperties.subDict("surfaceMasters");
    Info << "Reading surface masters \n" << endl;
    surfaceMastersNames = surfacesDict.toc();
}


speciesTable surfaceMasters
(
    surfaceMastersNames
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "Reading two-phase properties\n" << endl;
twoPhaseMixture twoPhaseProperties(U, phi, "alpha1");

Info<< "Building Interface instance\n" << endl;
interfacePropertiesAM interface(alpha1, U, twoPhaseProperties);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


Info<< "Reading field Surf\n" << endl;
volScalarField Surf
(
    IOobject
    (
        "Surf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "createConcMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    /*label patchID = mesh.boundaryMesh().findPatchID("Grainwalls");

    const polyPatch& cPatch = mesh.boundaryMesh()[patchID];

    const surfaceScalarField& magSf = mesh.magSf();

    scalar patchArea = 0.0;
    forAll(cPatch, facei)
    {
        patchArea += magSf.boundaryField()[patchID][facei];
    }*/

/*Info << "Reading initial component Succinate" << endl;
volScalarField C_s_init ("C_s_init", Y[0]);
forAll(C_s_init, i) {
		C_s_init[i] = alpha1[i]*scalar(0.2);	
	}
Y[0] = C_s_init;
Y[0].write();*/