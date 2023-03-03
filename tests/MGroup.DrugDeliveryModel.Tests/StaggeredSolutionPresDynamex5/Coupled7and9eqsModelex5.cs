using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.DrugDeliveryModel.Tests.EquationModels;
using MGroup.NumericalAnalyzers.Dynamic;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.Staggered;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.Constitutive.Structural;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.NumericalAnalyzers;
using MGroup.Solvers.Direct;
using Xunit;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Solution;
using MGroup.FEM.ConvectionDiffusion.Isoparametric;
using MGroup.FEM.Structural.Continuum;
using ab = MGroup.DrugDeliveryModel.Tests.Commons.Utilities;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
	public class Coupled7and9eqsModelex5
    {
        public Eq78ModelProviderForStaggeredSolutionex5 Eq78ModelProvider { get; set; }
        public Eq9ModelProviderForStaggeredSolutionex5 Eq9ModelProvider { get; set; }

        public Model[] model;



        //TODO put analysis time stepping where it belongs1 (pithanws sto Coupled7and9eqsSolution.cs h Coupled7and9eqsModel.cs)
        private GenericAnalyzerState[] analyzerStates, nlAnalyzerStates;
        private IParentAnalyzer[] parentAnalyzers;
        private IChildAnalyzer[] nlAnalyzers;
        private ISolver[] parentSolvers;


        public int CurrentTimeStep { get; set; }

        public GenericAnalyzerState[] AnalyzerStates => analyzerStates;
        public GenericAnalyzerState[] NLAnalyzerStates => nlAnalyzerStates;
        public IParentAnalyzer[] ParentAnalyzers => parentAnalyzers;
        public IChildAnalyzer[] NLAnalyzers => nlAnalyzers;
        public ISolver[] ParentSolvers => parentSolvers;
        public ComsolMeshReader Reader => reader;

        private ComsolMeshReader reader;

        private Dictionary<int, double> lambda;
        private Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints;
        private Dictionary<int, double[]> div_vs;

        private double timeStep;
        private double totalTime;

        private int incrementsPerStep;

        public Coupled7and9eqsModelex5(Eq78ModelProviderForStaggeredSolutionex5 eq78ModelProvider,
                                     Eq9ModelProviderForStaggeredSolutionex5 eq9ModelProvider, ComsolMeshReader comsolReader,
            Dictionary<int, double> lambda, Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints,
            Dictionary<int, double[]> div_vs, double timeStep, double totalTime, int incrementsPerStep)
        {
            Eq9ModelProvider = eq9ModelProvider;
            Eq78ModelProvider = eq78ModelProvider;
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;


            analyzerStates = new GenericAnalyzerState[2];
            nlAnalyzerStates = new GenericAnalyzerState[2];
            parentAnalyzers = new IParentAnalyzer[2];
            nlAnalyzers = new IChildAnalyzer[2];
            parentSolvers = new ISolver[2];

            reader = comsolReader;

            this.pressureTensorDivergenceAtElementGaussPoints = pressureTensorDivergenceAtElementGaussPoints;
            this.lambda = lambda;
            this.div_vs = div_vs;

            this.timeStep = timeStep;
            this.totalTime  = totalTime;
            this.incrementsPerStep = incrementsPerStep;

            // intialize array ofm models1.
            model = new Model[2];
        }


        public void CreateModel(IParentAnalyzer[] analyzers, ISolver[] solvers)
        {
            //---------------------------------------
            // WARNING: do not initialize shared dictionarys because they have been passed by refernce in ewuationModel bilders.
            //---------------------------------------


            // update Shared quantities of Coupled model
            //foreach (var elem in reader.ElementConnectivity)
            //{ 
            //    lambda[elem.Key]= lambda0;
            //}
            foreach (var elem in reader.ElementConnectivity)
            {
                //---------------
                //no coupling commented out update of shared quantities
                //
                pressureTensorDivergenceAtElementGaussPoints[elem.Key] = ab.ScalePressureTensorDiv(((ConvectionDiffusionElement3D)model[0].ElementsDictionary[elem.Key]).pressureTensorDivergenceAtGaussPoints);
            }
            foreach (var elem in reader.ElementConnectivity)
            {

                div_vs[elem.Key] = ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocityDivergence;
            }


            //create models with them 
            model = new Model[2];


            model[0] = Eq78ModelProvider.GetModel();
            Eq78ModelProvider.AddEq78ModelAppropriateBCs(model[0]);
            //Eq78ModelProvider.AddTopAndBottomBCsDistributedPeripheral(model[0]); // PROSOXH!!!! modify se duo shmeia tis BCs
            (analyzers[0], solvers[0], nlAnalyzers[0]) = Eq78ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //TODo
            model[1] = Eq9ModelProvider.GetModel();
            Eq9ModelProvider.AddBottomBCs(model[1]);
            //Eq9ModelProvider.AddEq9ModelAppropriateBCs(model[1]);
            Eq9ModelProvider.AddEq9ModelLoads(model[1]);
            (analyzers[1], solvers[1], nlAnalyzers[1]) = Eq9ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[1], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            for (int i = 0; i < analyzers.Length; i++)
            {
                analyzers[i].Initialize(true);
                if (analyzerStates[i] != null)
                {
                    analyzers[i].CurrentState = analyzerStates[i];
                }

                if (nlAnalyzerStates[i] != null)
                {
                    nlAnalyzers[i].CurrentState = nlAnalyzerStates[i];
                }
            }
        }

        public void CreateModelFirstTime(IParentAnalyzer[] analyzers, ISolver[] solvers)
        {
            //create models with initial 

            model = new Model[2];


            model[0] = Eq78ModelProvider.GetModel();
            Eq78ModelProvider.AddEq78ModelAppropriateBCs(model[0]);
            //Eq78ModelProvider.AddTopAndBottomBCsDistributedPeripheral(model[0]); // PROSOXH!!!! modify se duo shmeia tis BCs
            (analyzers[0], solvers[0], nlAnalyzers[0]) = Eq78ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //TODo
            model[1] = Eq9ModelProvider.GetModel();
            //Eq9ModelProvider.AddEq9ModelAppropriateBCs(model[1]);
            Eq9ModelProvider.AddBottomBCs(model[1]);
            Eq9ModelProvider.AddEq9ModelLoads(model[1]);
            (analyzers[1], solvers[1], nlAnalyzers[1]) = Eq9ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[1], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            for (int i = 0; i < analyzers.Length; i++)
            {
                analyzers[i].Initialize(true);
                if (analyzerStates[i] != null)
                {
                    analyzers[i].CurrentState = analyzerStates[i];
                }

                if (nlAnalyzerStates[i] != null)
                {
                    nlAnalyzers[i].CurrentState = nlAnalyzerStates[i];
                }
            }
        }




    }
}
