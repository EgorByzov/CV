using MathKernel;
using ObjectModel.LightDistributions;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;
using ViewModelCore.Light_Distributions;
using ViewModelCore.MaterialProperties;
using ViewModelCore.OpticalSceneEntities.OriginEntities;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{
    using RuleBook = RuleBook<LensAxisymTwoAsphericZAView>;
    using ValidationRule = DelegateValidationRule<LensAxisymTwoAsphericZAView>;
    using BindRule = PropertyBindRule<LensAxisymTwoAsphericZAView>;

    public class LensAxisymTwoAsphericZAView : SolidLensZAParamsView
    {
        private LensAxisymTwoAsphericZAParams modelObject;
        public override SolidLensZAParams ModelObject
        {
            get { return modelObject; }
        }
        private double rInn, rOut, skirtHeight;

        # region Properties
        public double RInn
        {
            get { return rInn; }
            set
            {
                SetProperty<LensAxisymTwoAsphericZAView, double>(ref rInn, value);
            }
        }
        public double ROut
        {
            get { return rOut; }
            set
            {
                SetProperty<LensAxisymTwoAsphericZAView, double>(ref rOut, value);
            }
        }
        public double SkirtHeight
        {
            get { return skirtHeight; }
            set
            {
                SetProperty<LensAxisymTwoAsphericZAView, double>(ref skirtHeight, value);
            }
        }
        # endregion

        //Constructors
        static LensAxisymTwoAsphericZAView()
        {
            InitRules();
        }
        internal LensAxisymTwoAsphericZAView(LensAxisymTwoAsphericZAParams ZA, MaterialView material)
        {
            modelObject = ZA;

            rInn = ZA.RInn;
            rOut = ZA.ROut;
            skirtHeight = ZA.SkirtHeight;
            Wavelength = ZA.Wavelength;

            SubscribeToModelObject<LensAxisymTwoAsphericZAParams, LensAxisymTwoAsphericZAView>(ZA);
        }
        public static LensAxisymTwoAsphericZAView DefaultZAView()
        {
            var material = ApplicationPreferences.MaterialsBase.DefaultMaterial();
            var defaultZA =
                new LensAxisymTwoAsphericZAParams(material.modelObject);
            return new LensAxisymTwoAsphericZAView(defaultZA, material);
        }
        //public
        public override OpticalSceneObject ComputeOpticalElement(OpticalEntityView source,
                                                         RequiredSimulationResultView reqSim)
        {
            var inten = ((SolidSourceView)source.OEOrigin).IntensityModel.IntensityDistribution;
            var cs = source.CS;
            var lens = new LensAxisymTwoAspheric(modelObject, (IntensityAxisym)inten.ModelObject,
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

            // bind
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
