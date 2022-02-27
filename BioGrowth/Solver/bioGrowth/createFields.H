
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

#include "createPhi.H"

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

Info << "Reading component A" << endl;
volScalarField A
(
    IOobject
    (
	"A",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading component B" << endl;
volScalarField B
(
    IOobject
    (
	"B",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading component C" << endl;
volScalarField C
(
    IOobject
    (
	"C",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading component bio" << endl;
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


dimensionedTensor DT
(
    transportProperties.lookup("DT")
);

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

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	//std::ostringstream oss;
	
/*	double* surfConcentration_1;
	double* surfConcentration_2;*/
	
	int ncells = mesh.cells().size();
	double* surfConcentration_2 = new double [ncells];
	double* surfArea = new double [ncells];

//	int nsurf = 2;
	
	surfArea           = (double *)malloc((size_t)(ncells * sizeof(double)));
//	surfConcentration_1  = (double *)malloc((size_t)(ncells * sizeof(double)));
//	surfConcentration_2  = (double *)malloc((size_t)(ncells * sizeof(double)));
	
	forAll(mesh.cells(),celli)
	{
		forAll(surfaceMasters, j)
		{
			double area = 0.0;
			//double mole = 0.0;
	        const volScalarField::GeometricBoundaryField& Surfbf = Surf.boundaryField();
			forAll(Surfbf,patchi)
			{
				if (Surfbf[patchi].type() == "reactingWall")
				{
					const reactingWallFvPatchScalarField& Surfcap = refCast<const reactingWallFvPatchScalarField>(Surfbf[patchi]);
					const labelList& cellOwner = Surfcap.patch().faceCells();
					const surfaceScalarField& magSf = mesh.magSf();
					const wordList& masters = Surfcap.get_surface_masters();
					const scalarList& density   = Surfcap.get_density();//mol/m^2
					forAll(masters,i)
					{
						if (masters[i]==surfaceMasters[j])
						{
							forAll(Surfbf[patchi],facei)
							{
								if (cellOwner[facei]==celli)
								{
									//mole+=density[i]*magSf.boundaryField()[patchi][facei] / mesh.V()[cellOwner[facei]];//mol/L
									area+=magSf.boundaryField()[patchi][facei] / mesh.V()[cellOwner[facei]];//m^2/m^3
								}
							}
						}
					}
				}
			}
			
			surfArea[celli]	= area;
			//Info << "surface area of cell " << celli <<"is: " << surfArea[celli] << endl;
		}
	}
	