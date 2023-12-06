using MathKernel;
using ObjectModel.LightDistributions;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;
using ViewModelCore.Light_Distributions;
using ViewModelCore.MaterialProperties;
using ViewModelCore.OpticalSceneEntities.OriginEntities;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{

    using RuleBook = RuleBook<LensFreeformSphDomeZAView>;
    using ValidationRule = DelegateValidationRule<LensFreeformSphDomeZAView>;
    using BindRule = PropertyBindRule<LensFreeformSphDomeZAView>;

    public class LensFreeformSphDomeZAView : SolidLensZAParamsView
    {
        //fields
        private double rInn, rOut, skirtHeight;
        private LensFreeformSphDomeZAParams modelObject;

        public override SolidLensZAParams ModelObject
        {
            get { return modelObject; }
        }

        # region Propertie
        public double RInn
        {
            get { return rInn; }
            set
            {
                SetProperty<LensFreeformSphDomeZAView, double>(ref rInn, value);
            }
        }
        public double ROut
        {
            get { return rOut; }
            set
            {
                SetProperty<LensFreeformSphDomeZAView, double>(ref rOut, value);
            }
        }
        public double SkirtHeight
        {
            get { return skirtHeight; }
            set
            {
                SetProperty<LensFreeformSphDomeZAView, double>(ref skirtHeight, value);
            }
        }
        # endregion

        static LensFreeformSphDomeZAView()
        {
            InitRules();
        }

        internal LensFreeformSphDomeZAView(LensFreeformSphDomeZAParams za, MaterialView material)
        {
            modelObject = za;
            rInn = za.RInn;
            rOut = za.ROut;
            skirtHeight = za.SkirtHeight;
            Wavelength = za.Wavelength;
            SubscribeToModelObject<LensFreeformSphDomeZAParams, LensFreeformSphDomeZAView>(za);
        }
        public static LensFreeformSphDomeZAView DefaultZAView()
        {
            var material = ApplicationPreferences.MaterialsBase.DefaultMaterial();
            var defaultZA =
                new LensFreeformSphDomeZAParams(material.modelObject);
            return new LensFreeformSphDomeZAView(defaultZA, material);
        }
        //public
        public override OpticalSceneObject ComputeOpticalElement(OpticalEntityView source,
                                                         RequiredSimulationResultView reqSim)
        {
            var inten = ((SolidSourceView)source.OEOrigin).IntensityModel.IntensityDistribution;
            var cs = source.CS;
            var lens = new LensFreeformSphDome(modelObject, (IntensityAxisym)inten.ModelObject,
                                                 reqSim.ModelObject, cs.ModelObject);
            return new OpticalSceneObject(lens);
        }

        //private
        private static void InitRules()
        {
            RuleBook.ValidationRules.Add(new ValidationRule(
                "RInn",
                "RInnOutOfRange",
                x => (Epsilon.CompareNumerics(x.RInn, 0) > 0
                    && Epsilon.CompareNumerics(x.RInn, x.ROut) < 0)));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "ROut",
                "ROutOutOfRange",
                x => (Epsilon.CompareNumerics(x.ROut, 0) > 0
                    && Epsilon.CompareNumerics(x.RInn, x.ROut) < 0)));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "SkirtHeight",
                "SkirtHeightOutOfRange",
                x => Epsilon.CompareNumerics(x.SkirtHeight, 0) > 0));
            RuleBook.BindRules.Add(new BindRule(
                "RInn",
                "RInn",
                x => x.modelObject.RInn = x.RInn,
                x => x.rInn = x.modelObject.RInn
                ));
            RuleBook.BindRules.Add(new BindRule(
                "ROut",
                "ROut",
                x => x.modelObject.ROut = x.ROut,
                x => x.rOut = x.modelObject.ROut
                ));
            RuleBook.BindRules.Add(new BindRule(
                "SkirtHeight",
                "SkirtHeight",
                x => x.modelObject.SkirtHeight = x.SkirtHeight,
                x => x.skirtHeight = x.modelObject.SkirtHeight
                ));

        }

    }
}
