#include "CONSTANTS.H"
#include "Vector.H"
#include "ParmParse.H"
#include "BoxIterator.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "AMRIO.H"

#include "UsingNamespace.H"

void fillData(Vector<LevelData<FArrayBox>* > a_data,
              const Vector<ProblemDomain>&   a_domains,
              const Vector<int>&             a_refRatios,
              const Vector<Real>&            a_dx,
              const RealVect&                a_domainOrigin,
              const int&                     a_finestLevel)
{
  CH_TIME("fillData");

  for (int lev = 0; lev <= a_finestLevel; lev++)
  {
    LevelData<FArrayBox>& levelData     = *(a_data[lev]);
    const DisjointBoxLayout& levelGrids = levelData.getBoxes();

    RealVect numCells = a_domains[lev].size();

    // data is cell-centered...
    RealVect offset = 0.5 * a_dx[lev] * RealVect::Unit + a_domainOrigin;

    DataIterator levelDit = levelGrids.dataIterator();
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      FArrayBox& thisData = levelData[levelDit];

      BoxIterator bit(thisData.box());
      for (bit.begin(); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect loc(iv);

        loc *= a_dx[lev];
        loc += offset;

        Real theta = loc[0];
        Real phi   = loc[1];
        Real psi   = loc[2];

        Real cth = cos(theta);
        Real sth = sin(theta);
        
        Real cph = cos(phi);
        Real sph = sin(phi);
        
        Real cps = cos(psi);
        Real sps = sin(psi);
        
        thisData(iv,0) = Abs( cps*cth+sps*sph*sth) + Abs(-cps*sth+sps*sph*cth);
        thisData(iv,1) = Abs(             cph*sth) + Abs(             cph*cth);
        thisData(iv,2) = Abs(-sps*cth+cps*sph*sth) + Abs( sps*sth+cps*sph*cth);
      }
    } // end loop over grids on this level
  } // end loop over levels
}

void setupGrids(Vector<DisjointBoxLayout>& a_grids,
                Vector<ProblemDomain>&     a_domains,
                Vector<int>&               a_refRatios,
                Vector<Real>&              a_dx,
                RealVect&                  a_domainOrigin,
                int&                       a_finestLevel)
{
  CH_TIME("setupGrids");

  a_finestLevel = 0;
  ParmParse ppGrids("grids");

  // get grid generation parameters
  int maxLevel;
  ppGrids.get("max_level",maxLevel);

  int maxBoxSize;
  ppGrids.get("max_box_size",maxBoxSize);

  int blockFactor;
  ppGrids.get("block_factor",blockFactor);

  Real fillRatio;
  ppGrids.get("fill_ratio",fillRatio);

  // note that there only need to be numLevels-1 refinement ratios
  a_refRatios.resize(maxLevel);
  ppGrids.getarr("ref_ratio",a_refRatios,0,maxLevel);

  bool isPeriodic[SpaceDim];
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    isPeriodic[dir] = false;
    // isPeriodic[dir] = true;
  }

  IntVect numCells;
  Vector<int> numCellsVect(SpaceDim);
  ppGrids.getarr("num_cells",numCellsVect,0,SpaceDim);
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    numCells[dir] = numCellsVect[dir];
  }

  Vector<Real> domainOriginVect(SpaceDim);
  ppGrids.getarr("domain_origin",domainOriginVect,0,SpaceDim);
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    a_domainOrigin[dir] = M_PI*domainOriginVect[dir];
  }

  RealVect domainSize = RealVect::Unit;
  Vector<Real> domainSizeVect(SpaceDim);
  ppGrids.getarr("domain_size",domainSizeVect,0,SpaceDim);
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    domainSize[dir] = M_PI*domainSizeVect[dir];
  }

  // resize dataholders
  int maxNumLevels = maxLevel + 1;

  a_grids.resize(maxNumLevels);
  a_domains.resize(maxNumLevels);
  a_dx.resize(maxNumLevels,-1);

  a_finestLevel = 0;

  a_dx[0] = domainSize[0]/numCells[0];

  IntVect domLo = IntVect::Zero;
  IntVect domHi = numCells - IntVect::Unit;

  ProblemDomain baseDomain(domLo,domHi,isPeriodic);
  a_domains[0] = baseDomain;

  // set up refined domains, etc
  for (int lev = 1; lev <= maxLevel; lev++)
    {
      a_domains[lev] = a_domains[lev-1];
      a_domains[lev].refine(a_refRatios[lev-1]);

      a_dx[lev] = a_dx[lev-1] / a_refRatios[lev-1];
    }

  Vector<Vector<Box> > vectBoxes(maxLevel+1);

  {
    CH_TIME("BaseGridCreation");

    // generate base level grids
    domainSplit(baseDomain,vectBoxes[0],maxBoxSize,blockFactor);

    pout() << "Number of boxes: " << vectBoxes[0].size() << endl;

    Vector<int> procAssign(vectBoxes[0].size(),0);
    LoadBalance(procAssign,vectBoxes[0]);

    DisjointBoxLayout baseGrids(vectBoxes[0],procAssign,baseDomain);
    a_grids[0] = baseGrids;
  }

  if (maxLevel > 0)
  {
    MayDay::Error("AMR not currently handled");
#if 0
    bool read_grids = false;
    ppGrids.query("read_in_grids",read_grids);

    // tag on grad(rhs)
    int bufferSize = 1;
    BRMeshRefine meshGen(a_amrDomains[0],
                         a_refRatios,
                         fillRatio,
                         blockFactor,
                         bufferSize,
                         maxBoxSize);

    // to be used by MeshRefine...
    Vector<Vector<Box> > oldMeshes(maxLevel+1);
    oldMeshes[0] = vectBoxes[0];
    for (int lev=1; lev<oldMeshes.size(); lev++)
    {
      oldMeshes[lev].push_back(a_amrDomains[lev].domainBox());
    }

    Real refineThresh;
    ppGrids.get("refine_threshold",refineThresh);

    Real threshSqr = refineThresh*refineThresh;

    bool moreLevels = true;
    while (moreLevels)
    {
      // tag based on grad(rhs)
      // first need to allocate RHS
      Vector<LevelData<FArrayBox>* > tempRHS(a_finestLevel+1,NULL);
      for (int lev=0; lev<= a_finestLevel; lev++)
        {
          // note that we add a ghost cell to simplify gradients
          tempRHS[lev] = new LevelData<FArrayBox>(a_amrGrids[lev],
                                                  1,IntVect::Unit);
        }

      setRHS(tempRHS,a_amrDomains,a_refRatios,a_amrDx,
             a_finestLevel);

      Vector<IntVectSet> tags(a_finestLevel+1);

      for (int lev=0; lev<a_finestLevel+1; lev++)
        {
          const DisjointBoxLayout& levelGrids = a_amrGrids[lev];
          const LevelData<FArrayBox>& levelRHS = *tempRHS[lev];
          IntVectSet& levelTags = tags[lev];

          // compute mag(gradient)
          DataIterator dit = levelGrids.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              const FArrayBox& rhsFab = levelRHS[dit];
              // local storage foer gradient
              FArrayBox gradFab(levelGrids[dit],1);
              gradFab.setVal(0.0);
              Real thisGrad;

              BoxIterator bit(levelGrids[dit]);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv=bit();
                  for (int dir=0; dir<SpaceDim; dir++)
                    {
                      // use mag(undivided gradient)
                      IntVect hi = iv + BASISV(dir);
                      IntVect lo = iv - BASISV(dir);
                      thisGrad = rhsFab(hi,0) - rhsFab(lo,0);
                      gradFab(iv,0) += (thisGrad*thisGrad);
                    } // end loop over directions
                } // end loop over cells

              //gradFab now has mag(grad*dx)^2

              // tag where mag(gradient) > tolerance^2
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  if (gradFab(iv,0) > threshSqr)
                    {
                      levelTags |= iv;
                    }
                } // end loop over cells
            } // end loop over grids on this level

        } // end loop over levels


      // call meshRefine.
      for (int lev=1; lev<=a_finestLevel; lev++)
        {
          oldMeshes[lev] = vectBoxes[lev];
        }

      int topLevel = a_finestLevel;
      int newFinestLevel =  meshGen.regrid(vectBoxes,
                                           tags,
                                           0,
                                           topLevel,
                                           oldMeshes);


      // define new grids if necessary and test to see if we're done
      if (newFinestLevel > a_finestLevel)
      {
        a_finestLevel = newFinestLevel;

        // setup new grid hierarchy
        for (int lev=1; lev<=a_finestLevel; lev++)
          {
            Vector<int> procAssign(vectBoxes[lev].size(),0);
            LoadBalance(procAssign,vectBoxes[lev]);
            DisjointBoxLayout levelGrids(vectBoxes[lev],
                                         procAssign,
                                         a_amrDomains[lev]);
            a_amrGrids[lev] = levelGrids;
          }
      }
      else
      {
        moreLevels = false;
      }

      if (a_finestLevel == maxLevel)
      {
        moreLevels = false;
      }

      // clean up before starting again
      for (int lev=0; lev<tempRHS.size(); lev++)
      {
        delete tempRHS[lev];
      }
    } // end while (moreLevels)
#endif
  }

  // fill in remaining levels with empty DisjointBoxLayouts
  for (int lev = a_finestLevel + 1; lev <= maxLevel; lev++)
  {
    a_grids[lev] = DisjointBoxLayout();
  }
}

void runProblem()
{
  // CH_TIME("runProblem");
  CH_TIME("runProblem");

  // set up grids
  Vector<DisjointBoxLayout> grids;
  Vector<ProblemDomain>     domains;
  Vector<int>               refRatios;
  Vector<Real>              dx;
  RealVect                  domainOrigin;
  int                       finestLevel;

  setupGrids(grids,domains,refRatios,dx,domainOrigin,finestLevel);

  // allocate, initialize data
  int numLevels = grids.size();
  Vector<LevelData<FArrayBox>* > data(numLevels,NULL);

  int nComps = 3;
  IntVect ghostCells = IntVect::Zero;

  for (int lev = 0; lev <= finestLevel; lev++)
  {
    CH_TIME("createLevels");
    const DisjointBoxLayout& levelGrids = grids[lev];
    data[lev] = new LevelData<FArrayBox>(levelGrids,nComps,ghostCells);
  }

#if 0
  fillData(data,domains,refRatios,dx,domainOrigin,finestLevel);

#ifdef CH_USE_HDF5
  numLevels = finestLevel + 1;
 
  string fname = "equalRows.";

  char suffix[30];
  sprintf(suffix,"%dd.hdf5",SpaceDim);
  fname += suffix;

  Vector<string> varNames(3);
  varNames[0] = "row1";
  varNames[1] = "row2";
  varNames[2] = "row3";

  Real bogusVal = 1.0;

  WriteAMRHierarchyHDF5(fname,
                        grids,
                        data,
                        varNames,
                        domains[0].domainBox(),
                        dx[0],
                        bogusVal,
                        bogusVal,
                        refRatios,
                        numLevels);
#endif // end if HDF5
#endif

  // clean up
  for (int lev = 0; lev < data.size(); lev++)
  {
    CH_TIME("deleteLevels");
    delete data[lev];
  }
}

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // scoping...
  {
    if (argc < 2)
      {
        cerr << endl;
        cerr << "Error: No input file specified..." << endl;
        cerr << endl;
        cerr << "Usage: " << argv[0] << " <input_file_name> " << endl;
        cerr << endl;

        exit(0);
      }

    char* in_file = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,in_file);

    runProblem();
  }
  //end scoping...
  
#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif

  return(0);
}
