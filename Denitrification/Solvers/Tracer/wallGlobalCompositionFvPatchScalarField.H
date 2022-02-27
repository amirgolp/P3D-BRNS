/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

#ifndef wallGlobalCompositionFvPatchScalarField_H
#define wallGlobalCompositionFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class wallGlobalCompositionFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class wallGlobalCompositionFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

protected:

    // Protected data

      word speciesName_;

      scalar H_;

public:

    //- Runtime type information
    TypeName("wallGlobalComposition");


    // Constructors

        //- Construct from patch and internal field
        wallGlobalCompositionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wallGlobalCompositionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  wallGlobalCompositionFvPatchScalarField onto a new patch
        wallGlobalCompositionFvPatchScalarField
        (
            const wallGlobalCompositionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wallGlobalCompositionFvPatchScalarField
        (
            const wallGlobalCompositionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new wallGlobalCompositionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        wallGlobalCompositionFvPatchScalarField
        (
            const wallGlobalCompositionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new wallGlobalCompositionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions



        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
