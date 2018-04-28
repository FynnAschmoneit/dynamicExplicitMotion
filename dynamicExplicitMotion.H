/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::dynamicExplicitMotion

Description
    

SourceFiles
    dynamicExplicitMotion.C

modified: Fynn Aschmoneit, February 2018, fyna@env.dtu.dk
\*---------------------------------------------------------------------------*/

#ifndef dynamicExplicitMotion_H
#define dynamicExplicitMotion_H

#include "dynamicFvMesh.H"
#include "dictionary.H"
#include "pointIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicExplicitMotion Declaration
\*---------------------------------------------------------------------------*/

class dynamicExplicitMotion
:
    public dynamicFvMesh
{
    // Private data

        dictionary dynamicMeshCoeffs_;

        scalar amplitude_;
        scalar frequency_;
        word  refPlaneDrc_;
        scalar refPlaneCrdZ_;
        scalar refPlaneCrdY_;
        word  refK1Drc_;
        word  refK2Drc_;
        scalar refK1Md_;
        scalar refK2Md_;

        pointIOField stationaryPoints_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicExplicitMotion(const dynamicExplicitMotion&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicExplicitMotion&);


public:

    //- Runtime type information
    TypeName("dynamicExplicitMotion");


    // Constructors

        //- Construct from IOobject
        dynamicExplicitMotion(const IOobject& io);


    //- Destructor
    ~dynamicExplicitMotion();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
