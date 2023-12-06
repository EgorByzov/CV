using System;
using System.IO;
using MathKernel;
using ObjectModel.Functions;
using ObjectModel.LightDistributions;
using ObjectModel.OpticalSceneEntities.OriginEntities;
using ViewModelCore.Common;
using ViewModelCore.Light_Distributions;
using ViewModelCore.Functions;
using ViewModelCore.IO;

namespace ViewModelCore.OpticalSceneEntities.OriginEntities
{
    using RuleBook = RuleBook<SolidSourceView>;
    using ValidationRule = DelegateValidationRule<SolidSourceView>;
    using BindRule = PropertyBindRule<SolidSourceView>;

    public class SolidSourceView : OpticalEntityOriginView, IEmittableSolidContainer, IIntensityModelContainer
    {
        # region Fields

        private SolidSource modelObject
        {
            get { return (SolidSource) base.ModelObject; }
        }

        private double flux;
        private int numRays;
        private EmittableSolidView solid;
        private IntensityModelView intensityModel;

        # endregion

        # region Properties

        #region IEmittableSolidContainer members

        public EmittableSolidType SolidType
        {
            get { return Solid.Type; }
            set
            {
                if (value == SolidType) return;

                Solid = EmittableSolidView.GetSolidByType(value);
            }
        }

        public EmittableSolidView Solid
        {
            get { return solid; }
            set
            {
                if (SetProperty(ref solid, value))
                {
                    modelObject.Solid = Solid.ModelObject;
                    OnPropertyChanged("SolidType");
                }
            }
        }

        #endregion

        #region IIntensityModelContainer members

        public IntensityModelType IntensityModelType
        {
            get { return IntensityModel.Type; }
            set
            {
                if (value == IntensityModelType) return;

                IntensityModel = IntensityModelView.GetModelByType(value);
            }
        }

        public IntensityModelView IntensityModel
        {
            get { return intensityModel; }
            set
            {
                if (SetProperty(ref intensityModel, value))
                {
                    modelObject.IntensityModel = IntensityModel.ModelObject;
                    OnPropertyChanged("IntensityModelType");
                }
            }
        }

        #endregion

        public double Flux
        {
            get { return flux; }
            set { SetProperty<SolidSourceView, double>(ref flux, value); }
        }

        public int NumRays
        {
            get { return numRays; }
            set { SetProperty<SolidSourceView, int>(ref numRays, value); }
        }
        public  DistributionProfileView[] Profiles
        {
            get
            {
                if(IntensityModel.ModelObject.IntensityDistribution is IntensityAxisym)
                {
                    var profiles = new DistributionProfileView[1];
                    profiles[0] = new DistributionProfileView(new IntensityDistributionProfile(0, IntensityModel.ModelObject.IntensityDistribution));
                    return profiles;}

                else
                {
                    var profiles = new DistributionProfileView[2];
                    profiles[0] = new DistributionProfileView(new IntensityDistributionProfile(0, IntensityModel.ModelObject.IntensityDistribution));
                    profiles[1] = new DistributionProfileView(new IntensityDistributionProfile(Math.PI/2, IntensityModel.ModelObject.IntensityDistribution));
                    return profiles;
                }
            }}
        # endregion

        public SolidSourceView() : this(new SolidSource())
        {
        }

        public SolidSourceView(SolidSource solidSource)
            : base(solidSource)
        {
            flux = solidSource.Flux;
            numRays = solidSource.NumRays;

            Solid = EmittableSolidView.Wrap(modelObject.Solid);
            IntensityModel = IntensityModelView.Wrap(modelObject.IntensityModel);

            SubscribeToModelObject<SolidSource, SolidSourceView>(solidSource);
        }

        static SolidSourceView()
        {
            InitRules();
        }

        private static void InitRules()
        {
            RuleBook.ValidationRules.Add(new ValidationRule(
                "Flux",
                "FluxOutOfRange",
                x => Epsilon.CompareNumerics(x.Flux, 0) > 0)
                );
            RuleBook.ValidationRules.Add(new ValidationRule(
                "NumRays",
                "NumRaysOutOfRange",
                x => x.NumRays > 0)
                );
            RuleBook.BindRules.Add(new BindRule
                (
                "Flux",
                "Flux",
                x => x.modelObject.Flux = x.Flux,
                x => x.Flux = x.modelObject.Flux
                ));
            RuleBook.BindRules.Add(new BindRule
                (
                "NumRays",
                "NumRays",
                x => x.modelObject.NumRays = x.NumRays,
                x => x.NumRays = x.modelObject.NumRays
                ));
        }



        public void SetLambertIntensityModel( )
        {
            var intenAx = new IntensityAxisym(new Cosinus(), new SphericalRing(0, Math.PI/2));
            var intenDistr = new ContiniousIntensityModel(intenAx);
            IntensityModel = new ContiniousIntensityModelView(intenDistr);
        }

        public void SetIsotropicIntensityModel()
        {
            var intenAx = new IntensityAxisym(new Constant(1), new SphericalRing(0, Math.PI / 2));
            var intenDistr = new ContiniousIntensityModel(intenAx);
            IntensityModel = new ContiniousIntensityModelView(intenDistr);
        }

        public void SetGaussIntensityModel()
        {
            var intenAx = new IntensityAxisym(new Gauss2D(), new SphericalRing(0, Math.PI / 2));
            var intenDistr = new ContiniousIntensityModel(intenAx);
            IntensityModel = new ContiniousIntensityModelView(intenDistr);
        }

        public void SetUserIntensityModel(Functions.ProfileView profile)
        {
            var intenAx = new IntensityAxisym(profile.Func2D.ModelObject, new SphericalRing(0, profile.MaxArg));
            var intenDistr = new ContiniousIntensityModel(intenAx);
            IntensityModel = new ContiniousIntensityModelView(intenDistr);
        }

        public void SaveIntenDistr(string fileName, LightDistributionFileType type, IESOptionsView options)
        {
            var distr = IntensityModel.IntensityDistribution;
            distr.SaveToFile(fileName, type, options);
        }

        public void OpenIntenDistrFromRayfile(string fileName)
        {
            IntensityModel = new DiscreteIntensityModelView(new DiscreteIntensityModel());
            ((DiscreteIntensityModelView) IntensityModel).FileName = fileName;
        }

        public void OpenIntenDistrFromIES(string fileName)
        {
            var distr = InputOutputView.LoadIntensityDistributionFromIES(fileName);
            IntensityModel = new ContiniousIntensityModelView(distr);
        }

        public void OpenIntenDistrFromInnerFormatFile(string fileName)
        {
        }

        public Functions.ProfileView TransformIntenDistrProfile2Function()
        {
            Functions.ProfileView profile;
            if (IntensityModel.ModelObject.IntensityDistribution is IntensityAxisym)
            {
                var func = ((IntensityAxisym) IntensityModel.ModelObject.IntensityDistribution).Profile;
                profile = new Functions.ProfileView(func);
            }
            //TODO
            else //(IntensityModel is DiscreteIntensityModelView)
            {
                profile = IntensityModel.IntensityDistribution.Profiles[0].Convert2ProfileView(0, Math.PI/2);
            }
            profile.Func2D = profile.Func2D.GetEditableFunc(profile.MinArg, profile.MaxArg);
            return profile;
        }

        public string GetTypeIntensityModel()
        {
            var model = IntensityModel.ModelObject;
            if (model is DiscreteIntensityModel)
            {
                var discrete = (DiscreteIntensityModelView) IntensityModel;
                var extension = Path.GetExtension(discrete.FileName);
                if (extension == null) return "lblFile";
                if (extension.Contains("ray") || extension.Contains("dat"))
                    return "lblRayfile";

                return "lblFile";
            }
            else
            {
                if (model.IntensityDistribution is IntensityByFunc3D) return "lblIESfile";
                var func = ((IntensityAxisym) model.IntensityDistribution).Profile.Func;
                if (func is Cosinus)
                    return "lblLambert";
                if (func is Constant)
                    return "lblIsotropic";
                if (func is Gauss2D)
                    return "lblGauss";

                return "lblUserDefined";
            }
        }

        public bool IsAxisymIntensityModel()
        {
            return IntensityModel.ModelObject.IntensityDistribution is IntensityAxisym;
        }
    }
}
