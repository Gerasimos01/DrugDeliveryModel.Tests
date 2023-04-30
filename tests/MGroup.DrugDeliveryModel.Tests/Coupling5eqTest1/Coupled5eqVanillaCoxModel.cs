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
using MGroup.DrugDeliveryModel.Tests.PreliminaryModels;
using MGroup.Solvers.AlgebraicModel;
using MGroup.LinearAlgebra.Matrices;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class Coupled5eqVanillaCoxModel
    {
        public Eq78ModelProviderForStaggeredSolutionex7ref Eq78ModelProvider { get; set; }
        public CoxVanillaSourceModelBuilder CoxModelProvider { get; set; }
        public Eq9ModelProviderForStaggeredSolutionEx7Ref Eq9ModelProvider { get; set; }
        public TCellModelProvider TCellModelProvider { get; }

        public DistributedOdeModelBuilder DistributedOdeModelProvider { get; set; }

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
        private Dictionary<int, double[]> FluidSpeed;
        private double kth;//TODO5eq mpws uparxoun duo kth?
        private Dictionary<int, double[][]> SolidVelocityAtElementGaussPoints;
        private Dictionary<int, double> DomainCOx { get; }

        private readonly Dictionary<int, double> domainT; // [cells]

        private double timeStep;
        private double totalTime;

        private int incrementsPerStep;

        public Coupled5eqVanillaCoxModel(Eq78ModelProviderForStaggeredSolutionex7ref eq78ModelProvider, CoxVanillaSourceModelBuilder coxModelProvider,
                                     Eq9ModelProviderForStaggeredSolutionEx7Ref eq9ModelProvider, TCellModelProvider tCellModelProvider,
                                     DistributedOdeModelBuilder distributedOdeModelProvider, ComsolMeshReader comsolReader,
                                      Dictionary<int, double> domainCOx, Dictionary<int, double> T,
            Dictionary<int, double> lambda, Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints,
            Dictionary<int, double[]> div_vs, Dictionary<int, double[]> FluidSpeed, Dictionary<int, double[][]> solidVelocityAtElementGaussPoints,
            double kth, double timeStep, double totalTime, int incrementsPerStep)
        {
            Eq9ModelProvider = eq9ModelProvider;
            CoxModelProvider = coxModelProvider;
            Eq78ModelProvider = eq78ModelProvider;
            TCellModelProvider = tCellModelProvider;
            DistributedOdeModelProvider= distributedOdeModelProvider;

            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;


            analyzerStates = new GenericAnalyzerState[5];
            nlAnalyzerStates = new GenericAnalyzerState[5];
            parentAnalyzers = new IParentAnalyzer[5];
            nlAnalyzers = new IChildAnalyzer[5];
            parentSolvers = new ISolver[5];

            reader = comsolReader;

            this.pressureTensorDivergenceAtElementGaussPoints = pressureTensorDivergenceAtElementGaussPoints;
            this.lambda = lambda;
            this.div_vs = div_vs;
            this.FluidSpeed = FluidSpeed;
            this.kth = kth;
            this.DomainCOx = domainCOx;
            this.domainT = T;


            this.SolidVelocityAtElementGaussPoints = solidVelocityAtElementGaussPoints;

            this.timeStep = timeStep;
            this.totalTime = totalTime;
            this.incrementsPerStep = incrementsPerStep;

            // intialize array ofm models1.
            model = new Model[5];
        }


        public void CreateModel(IParentAnalyzer[] analyzers, ISolver[] solvers)
        {
            //---------------------------------------
            // WARNING: do not initialize shared dictionarys because they have been passed by refernce in ewuationModel bilders.
            //---------------------------------------


            // update Shared quantities of Coupled model
            DistributedOdeModelProvider.RetrieveLambdaSolution(ParentSolvers[4], NLAnalyzers[4], model[4], lambda);
            CoxModelProvider.UpdateGausspointValuesOfElements(DomainCOx, ParentSolvers[2], NLAnalyzers[2], model[2], CoxModelProvider.algebraicModel);
            TCellModelProvider.UpdateGausspointValuesOfElements(domainT, ParentSolvers[3], NLAnalyzers[3], model[3], TCellModelProvider.algebraicModel);
            //foreach (var elem in reader.ElementConnectivity)
            //{ 
            //    lambda[elem.Key]= lambda0;
            //}
            foreach (var elem in reader.ElementConnectivity)
            {
                pressureTensorDivergenceAtElementGaussPoints[elem.Key] = ((ConvectionDiffusionElement3D)model[0].ElementsDictionary[elem.Key]).pressureTensorDivergenceAtGaussPoints;
            }
            foreach (var elem in reader.ElementConnectivity)
            {
                div_vs[elem.Key] = ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocityDivergence;
            }
            foreach (var elem in reader.ElementConnectivity)
            {//CALCULATE vf = kP + vs
                FluidSpeed[elem.Key] = new double[3]
                {
                    (-(pressureTensorDivergenceAtElementGaussPoints[elem.Key][0][0] * kth) + ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocity[0][0]) * 1000,
                    (-(pressureTensorDivergenceAtElementGaussPoints[elem.Key][0][1] * kth) + ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocity[0][1]) * 1000,
                    (-(pressureTensorDivergenceAtElementGaussPoints[elem.Key][0][2] * kth) + ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocity[0][2]) * 1000,
                };
            }
            foreach (var elem in reader.ElementConnectivity)
            {
                var velocityAtGP0 = ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocity[0];
                SolidVelocityAtElementGaussPoints[elem.Key][0][0] = velocityAtGP0[0] * 1000;
                SolidVelocityAtElementGaussPoints[elem.Key][0][1] = velocityAtGP0[1] * 1000;
                SolidVelocityAtElementGaussPoints[elem.Key][0][2] = velocityAtGP0[2] * 1000;
            }


            model = new Model[5];

            //Create model for eq78 (fluid pressure)
            model[0] = Eq78ModelProvider.GetModel();
            Eq78ModelProvider.AddBoundaryConditions(model[0]);
            (analyzers[0], solvers[0], nlAnalyzers[0]) = Eq78ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //Create model for eq9 (hyper-elastic material)
            model[1] = Eq9ModelProvider.GetModel();
            Eq9ModelProvider.AddBoundaryConditions(model[1]);
            (analyzers[1], solvers[1], nlAnalyzers[1]) = Eq9ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[1], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //Create model for eq13 (cox)
            model[2] = CoxModelProvider.GetModel();
            CoxModelProvider.AddBoundaryConditions(model[2]);
            (analyzers[2], solvers[2], nlAnalyzers[2]) = CoxModelProvider.GetAppropriateSolverAnalyzerAndLog(model[2], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //Create model for T equation 
            model[3] = TCellModelProvider.GetModel();
            TCellModelProvider.AddBoundaryConditions(model[3]);
            (analyzers[3], solvers[3], nlAnalyzers[3]) = TCellModelProvider.GetAppropriateSolverAnalyzerAndLog(model[3], timeStep, totalTime, CurrentTimeStep);

            //Create model for lambda equation 
            model[4] = DistributedOdeModelProvider.GetModel();
            DistributedOdeModelProvider.AddBoundaryConditions(model[4]);
            if (CurrentTimeStep == 0)
            {
                DistributedOdeModelProvider.AddInitialConditionsForTheRestOfBulkNodes(model[4]);

            }
            (analyzers[4], solvers[4], nlAnalyzers[4]) = DistributedOdeModelProvider.GetAppropriateSolverAnalyzerAndLog(model[4], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

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
            if (!(CurrentTimeStep == 0))
            {
                DistributedOdeModelProvider.RetrieveLambdaSolution(ParentSolvers[4], NLAnalyzers[4], model[4], lambda);
                CoxModelProvider.UpdateGausspointValuesOfElements(DomainCOx, ParentSolvers[2], NLAnalyzers[2], model[2], CoxModelProvider.algebraicModel);
                TCellModelProvider.UpdateGausspointValuesOfElements(domainT, ParentSolvers[3], NLAnalyzers[3], model[3], TCellModelProvider.algebraicModel);

                foreach (var elem in reader.ElementConnectivity)
                {
                    pressureTensorDivergenceAtElementGaussPoints[elem.Key] = ((ConvectionDiffusionElement3D)model[0].ElementsDictionary[elem.Key]).pressureTensorDivergenceAtGaussPoints;
                }
                foreach (var elem in reader.ElementConnectivity)
                {
                    div_vs[elem.Key] = ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocityDivergence;
                }
                foreach (var elem in reader.ElementConnectivity)
                {//CALCULATE vf = kP + vs
                 //vs
                    FluidSpeed[elem.Key] = new double[3]
                    {
                         (-(pressureTensorDivergenceAtElementGaussPoints[elem.Key][0][0] * kth) + ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocity[0][0]) * 1000,
                         (-(pressureTensorDivergenceAtElementGaussPoints[elem.Key][0][1] * kth) + ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocity[0][1]) * 1000,
                         (-(pressureTensorDivergenceAtElementGaussPoints[elem.Key][0][2] * kth) + ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocity[0][2]) * 1000,
                    };
                }
                foreach (var elem in reader.ElementConnectivity)
                {
                    var velocityAtGP0 = ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocity[0];
                    SolidVelocityAtElementGaussPoints[elem.Key][0][0] = velocityAtGP0[0] * 1000;
                    SolidVelocityAtElementGaussPoints[elem.Key][0][1] = velocityAtGP0[1] * 1000;
                    SolidVelocityAtElementGaussPoints[elem.Key][0][2] = velocityAtGP0[2] * 1000;
                }

            }


            model = new Model[5];

            //Create Initial Model eq78 (fluid pressure)
            model[0] = Eq78ModelProvider.GetModel();
            Eq78ModelProvider.AddBoundaryConditions(model[0]);
            if (CurrentTimeStep == 0)
            {
                //Eq78ModelProvider.AddEq78ModelInitialConditions(model[0]);
            }
            (analyzers[0], solvers[0], nlAnalyzers[0]) = Eq78ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //Create model for eq9 (hyperelastic material)
            model[1] = Eq9ModelProvider.GetModel();
            Eq9ModelProvider.AddBoundaryConditions(model[1]);
            (analyzers[1], solvers[1], nlAnalyzers[1]) = Eq9ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[1], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //Create model for eq8 (Cox)
            model[2] = CoxModelProvider.GetModel();
            CoxModelProvider.AddBoundaryConditions(model[2]);
            (analyzers[2], solvers[2], nlAnalyzers[2]) = CoxModelProvider.GetAppropriateSolverAnalyzerAndLog(model[2], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //Create model for T equation (cox)
            model[3] = TCellModelProvider.GetModel();
            TCellModelProvider.AddBoundaryConditions(model[3]);
            (analyzers[3], solvers[3], nlAnalyzers[3]) = TCellModelProvider.GetAppropriateSolverAnalyzerAndLog(model[3], timeStep, totalTime, CurrentTimeStep);

            //Create model for lambda equation 
            model[4] = DistributedOdeModelProvider.GetModel();
            DistributedOdeModelProvider.AddBoundaryConditions(model[4]);
            if (CurrentTimeStep == 0)
            {
                DistributedOdeModelProvider.AddInitialConditionsForTheRestOfBulkNodes(model[4]);

            }
            (analyzers[4], solvers[4], nlAnalyzers[4]) = DistributedOdeModelProvider.GetAppropriateSolverAnalyzerAndLog(model[4], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);


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

        public void SaveStateFromElements()
        {
            Eq9ModelProvider.SaveStateFromElements(model[1]);
        }

        
    }
}
