// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Info << "Reading solution species dictionary" << endl;
const dictionary& solutionSpeciesDict = twoPhaseProperties.subDict("solutionSpecies");
Info << "Reading solution species name\n" << endl;
speciesTable solutionSpecies
(
	solutionSpeciesDict.toc()
);

Info << "Create composition vector\n" << endl;
PtrList<volScalarField> Y(solutionSpecies.size());

forAll(solutionSpecies, i)
{
	IOobject header
	(
		solutionSpecies[i],
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ
	);

	// check if field exists and can be read
	if (header.headerOk())
	{
		Y.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					solutionSpecies[i],
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			)
		);
	}
	else
	{
		volScalarField Ydefault
		(
			IOobject
			(
				"Ydefault",
				mesh.time().timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			),
			mesh
		);

		Y.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					solutionSpecies[i],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				Ydefault
			)
		);
	}
	mesh.schemesDict().setFluxRequired(solutionSpecies[i]);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info << "Read species diffusion coefficients\n" << endl;
PtrList<dimensionedScalar> D1Y(solutionSpecies.size());
PtrList<dimensionedScalar> D2Y(solutionSpecies.size());
forAll(solutionSpecies, i)
{
	const dictionary& subdict = solutionSpeciesDict.subDict(solutionSpecies[i]);
	D1Y.set
	(
		i,
		new dimensionedScalar(subdict.lookup("D1"))
	);
	D2Y.set
	(
		i,
		new dimensionedScalar(subdict.lookup("D2"))
	);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info << "Read species molecular weight\n" << endl;
PtrList<dimensionedScalar> Mw(solutionSpecies.size());

forAll(solutionSpecies, i)
{
	const dictionary& subdict = solutionSpeciesDict.subDict(solutionSpecies[i]);
	Mw.set
	(
		i,
		new dimensionedScalar(subdict.lookup("Mw"))
	);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info << "Read species Henry's constant\n" << endl;
PtrList<dimensionedScalar> HY(solutionSpecies.size());

forAll(solutionSpecies, i)
{
	const dictionary& subdict = solutionSpeciesDict.subDict(solutionSpecies[i]);
	HY.set
	(
		i,
		new dimensionedScalar(subdict.lookup("H"))
	);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Info << "Create composition vector in each phase\n" << endl;
/*PtrList<volScalarField> Y1(solutionSpecies.size());
PtrList<volScalarField> Y2(solutionSpecies.size());

forAll(solutionSpecies, i)
{
	volScalarField& Yi = Y[i];
	dimensionedScalar HYi = HY[i];
	Y1.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"Y1_" + solutionSpecies[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			Yi/(alpha1+HYi*(1-alpha1))
		)
	);
	Y2.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"Y2_" + solutionSpecies[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			HYi*Yi/(alpha1+HYi*(1-alpha1))
		)
	);
}
*/
//**************************************************************************
scalarList maxY(solutionSpecies.size());
forAll(solutionSpecies,i)
{
	maxY[i] = readScalar(mesh.solutionDict().subDict("PISO").subDict("maxConcentration").lookup(solutionSpecies[i]));
}

Info << "Reading compressive flux coefficient\n" << endl;	
scalar cYi(readScalar(mesh.solutionDict().subDict("PISO").lookup("cYi")));

/*Info << "Reading equilibrium constant for SOM-DOM dynamics\n" << endl;	
dimensionedScalar kD_SOMDOM
(
    transportProperties.lookup("kD_SOMDOM")
);*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
Info << "Reading microbial species dictionary" << endl;
const dictionary& microbialSpeciesDict = twoPhaseProperties.subDict("microbialSpecies");
Info << "Reading microbial species name\n" << endl;
speciesTable microbialSpecies
(
	microbialSpeciesDict.toc()
);

Info << "Create composition vector\n" << endl;
PtrList<volScalarField> M(microbialSpecies.size());

forAll(microbialSpecies, i)
{
	IOobject header
	(
		microbialSpecies[i],
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ
	);

	// check if field exists and can be read
	if (header.headerOk())
	{
		Y.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					microbialSpecies[i],
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			)
		);
	}
	else
	{
		volScalarField Ydefault
		(
			IOobject
			(
				"Ydefault",
				mesh.time().timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			),
			mesh
		);

		Y.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					microbialSpecies[i],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				Ydefault
			)
		);
	}
	mesh.schemesDict().setFluxRequired(microbialSpecies[i]);
}
*/