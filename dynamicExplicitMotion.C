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

\*---------------------------------------------------------------------------*/

#include "dynamicExplicitMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicExplicitMotion, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicExplicitMotion, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicExplicitMotion::dynamicExplicitMotion(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    ),
    amplitude_(readScalar(dynamicMeshCoeffs_.lookup("amplitude"))),
    frequency_(readScalar(dynamicMeshCoeffs_.lookup("frequency"))),
    refPlaneDrc_(word(dynamicMeshCoeffs_.lookup("refPlaneDrc"))),
    refPlaneCrdZ_(readScalar(dynamicMeshCoeffs_.lookup("refPlaneCrdZ"))),
    refPlaneCrdY_(readScalar(dynamicMeshCoeffs_.lookup("refPlaneCrdY"))),
    refK1Drc_(word(dynamicMeshCoeffs_.lookup("refK1Drc"))),
    refK2Drc_(word(dynamicMeshCoeffs_.lookup("refK2Drc"))),
    refK1Md_(readScalar(dynamicMeshCoeffs_.lookup("refK1Md"))),
    refK2Md_(readScalar(dynamicMeshCoeffs_.lookup("refK2Md"))),
    stationaryPoints_
    (
        IOobject
        (
            "points",
            io.time().constant(),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{
    Info<< "Performing a dynamic mesh calculation: " << endl
        << "amplitude: " << amplitude_
        << " frequency: " << frequency_
        << " refPlaneDrc: " << refPlaneDrc_
        << " refPlaneCrdZ: " << refPlaneCrdZ_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicExplicitMotion::~dynamicExplicitMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicExplicitMotion::update()
{
    pointField newPoints = stationaryPoints_;
    vector::components displacementDirection;
    if (refPlaneDrc_ == "X" || refPlaneDrc_ == "x" ) displacementDirection = vector::X ;
    else if (refPlaneDrc_ == "Y" || refPlaneDrc_ == "y" ) displacementDirection = vector::Y ;
    else if (refPlaneDrc_ == "Z" || refPlaneDrc_ == "z" ) displacementDirection = vector::Z ;
    else {
	 displacementDirection = vector::X ;
	 Info<< "No reference plane orientation defined. Set refPlaneDrc to [X,x,Y,y,Z,z]. Default driection X applied." << endl;
    }
    
    // automatize like above for displacementDirection
    vector::components prpDrc1(vector::X);
    vector::components prpDrc2(vector::Y);
    
    Foam::Field<double> statPtDsplCmp(stationaryPoints_.component(displacementDirection));
    Foam::Field<double> statPtMd1Cmp(stationaryPoints_.component(prpDrc1));
    Foam::Field<double> statPtMd2Cmp(stationaryPoints_.component(prpDrc2));

    scalar waveLength1 = mag(max(statPtMd1Cmp)-min(statPtMd1Cmp));
    scalar waveLength2 = mag(mag(refPlaneCrdY_)-min(statPtMd2Cmp));

    Foam::Field<double> dsplRstr( (1.0-mag(statPtDsplCmp/max(statPtDsplCmp)))*(1.0-(statPtMd2Cmp-min(statPtMd2Cmp))/waveLength2) );   

    scalar ampl = amplitude_*sin( constant::mathematical::twoPi*frequency_*time().value() );

     newPoints.replace
    (
        displacementDirection,
	statPtDsplCmp			// current crd in displ direction
	+dsplRstr                       // field of maximal displacement
	*ampl				// time dependent ampitude
	*sin(constant::mathematical::twoPi*refK1Md_/waveLength1*statPtMd1Cmp) // sinosoidal displacement
    );

   fvMesh::movePoints(newPoints);

    volVectorField& U =
        const_cast<volVectorField&>(lookupObject<volVectorField>("U"));
    U.correctBoundaryConditions();

    return true;
}


// ************************************************************************* //
