
Info << "Cutting Plane" << endl;

// Get point and Plane
point xyz(0,2,0);
vector dir(0,1,0);

plane pl1(xyz,dir);

cuttingPlane cutPlane(pl1,mesh);

const labelList& cutCells = cutPlane.cutCells();
const word& name("someCells");

//scalarIOList z_plane
//(
//    IOobject
//    (
//        "z_plane",
//        runTime().timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::NO_WRITE
//    ),
//    scalarListList(cutCells.size())
//);


//forAll(cutCells,cellI)
//{
//   Info << "Z: " << z_plane[cellI] << endl;
//   Info << "C: " << mesh.C()[cutCells[cellI]] << endl;
   //z_plane[cellI]=mesh.C()[cutCells[cellI]];
//}

//z_plane.write();


//labelHashSet cCells(cutCells.size());
labelHashSet cCells;

forAll(cutCells,cellI)
{
        cCells.insert(cutCells[cellI]);
}

//volVectorField myField = cutPlane.sample(U_);

cellSet someCells(mesh,name,cCells);
//someCells.write();


fileName outputFile("p.csv");
OFstream os(db().time().timePath()/outputFile);
os << "First line of file.\n";
os << cutCells.size();
os << endl;
forAll(cutCells,i)          
{
os << mesh.C()[cutCells[i]];
os << p_;
os << endl;
}

//fvMeshSubset plane
//(
//    IOobject
//    (
//        "plane",
//        runTime().timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh
//);


//label patchI = -1;
//plane.setLargeCellSubset(someCells,patchI,true);
//
//
//wordList scalarNames(1);
//scalarNames[0] = "p";
//PtrList<volScalarField> scalarFlds(scalarNames.size()); 
//
//scalarFlds.set(0, plane.interpolate(p_));
//
//plane.subMesh().write();
//scalarFlds[0].write();
