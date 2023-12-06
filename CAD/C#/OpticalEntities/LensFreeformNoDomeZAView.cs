using MathKernel;
using ObjectModel.LightDistributions;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;
using ViewModelCore.Light_Distributions;
using ViewModelCore.MaterialProperties;
using ViewModelCore.OpticalSceneEntities.OriginEntities;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{
    using RuleBook = RuleBook<LensFreeformNoDomeZAView>;
    using ValidationRule = DelegateValidationRule<LensFreeformNoDomeZAView>;
    using BindRule = PropertyBindRule<LensFreeformNoDomeZAView>;

    public class LensFreeformNoDomeZAView : SolidLensZAParamsView
    {
        private LensFreeformNoDomeZAParams modelObject;

        public override SolidLensZAParams ModelObject
        {
            get { return modelObject; }
        }
        private double r0, skirtHeight;

        // Properties
        public double R0
        {
            get { return r0; }
            set
            {
                SetProperty<LensFreeformNoDomeZAView, double>(ref r0, value);
            }
        }

        public double SkirtHeight
        {
            get { return skirtHeight; }
            set
            {
                SetProperty<LensFreeformNoDomeZAView, double>(ref skirtHeight, value);
            }
        }

        //Constructors
        static LensFreeformNoDomeZAView()
        {
            InitRules();
        }

        internal LensFreeformNoDomeZAView(LensFreeformNoDomeZAParams ZA, MaterialView material)
        {
            modelObject = ZA;
            r0 = ZA.R0;
            skirtHeight = ZA.SkirtHeight;
            Wavelength = ZA.Wavelength;
            Material = material;
            SubscribeToModelObject<LensFreeformNoDomeZAParams, LensFreeformNoDomeZAView>(ZA);
        }

        public static LensFreeformNoDomeZAView DefaultZAView()
        {
            var material = ApplicationPreferences.MaterialsBase.DefaultMaterial();
            var defaultZA =
                new LensFreeformNoDomeZAParams(material.modelObject);
            return new LensFreeformNoDomeZAView(defaultZA, material);
        }

        //public
        public override OpticalSceneObject ComputeOpticalElement(OpticalEntityView source,
                                                         RequiredSimulationResultView reqSim)
        {
            var inten = ((SolidSourceView)source.OEOrigin).IntensityModel.IntensityDistribution;
            var cs = source.CS;
            var lens = new LensFreeformNoDome(modelObject, (IntensityAxisym)inten.ModelObject,
                                                 reqSim.ModelObject, cs.ModelObject);
            return new OpticalSceneObject(lens);
        }

        private static void InitRules()
        {
            RuleBook.ValidationRules.Add(new ValidationRule(
                "R0", "R0OutOfRange", x => Epsilon.CompareNumerics(x.R0, 0) > 0 ));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "SkirtHeight", "SkirtHeightOutOfRange", x => Epsilon.CompareNumerics(x.SkirtHeight, 0) > 0));
            RuleBook.BindRules.Add(new BindRule(
                 "R0",
                 "R0",
                 x => x.modelObject.R0 = x.R0,
                 x => x.r0 = x.modelObject.R0
                 ));
            RuleBook.BindRules.Add(new BindRule(
                "SkirtHeight",
                "SkirtHeight",
                x => x.modelObject.SkirtHeight = x.SkirtHeight,
                x => x.skirtHeight = x.modelObject.SkirtHeight
                ));
            RuleBook.BindRules.Add(new BindRule(
                "Wavelength",
                "Wavelength",
                x => x.modelObject.Wavelength = x.Wavelength,
                x => x.Wavelength = x.modelObject.Wavelength
                ));
        }

    }
}
