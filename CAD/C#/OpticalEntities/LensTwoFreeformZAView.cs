using MathKernel;
using ObjectModel.LightDistributions;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;
using ViewModelCore.Light_Distributions;
using ViewModelCore.MaterialProperties;
using ViewModelCore.OpticalSceneEntities.OriginEntities;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{
    using RuleBook = RuleBook<LensTwoFreeformZAView>;
    using ValidationRule = DelegateValidationRule<LensTwoFreeformZAView>;
    using BindRule = PropertyBindRule<LensTwoFreeformZAView>;

    public class LensTwoFreeformZAView : SolidLensZAParamsView
    {
        #region Fields

        private LensTwoFreeformZAParams modelObject;

        public override SolidLensZAParams ModelObject
        {
            get { return modelObject; }
        }
        private double rInn, rOut, skirtHeight, gamma0, hardOut;
        private MaterialView material;
        #endregion

        #region Properties
        public double RInn
        {
            get { return rInn; }
            set { SetProperty<LensTwoFreeformZAView, double>(ref rInn, value); }
        }

        public double ROut
        {
            get { return rOut; }
            set { SetProperty<LensTwoFreeformZAView, double>(ref rOut, value); }
        }

        public double SkirtHeight
        {
            get { return skirtHeight; }
            set { SetProperty<LensTwoFreeformZAView, double>(ref skirtHeight, value); }
        }

        public double Gamma0
        {
            get { return gamma0; }
            set { SetProperty<LensTwoFreeformZAView, double>(ref gamma0, value); }
        }

        public double HardOut
        {
            get {return hardOut; }
            set { SetProperty<LensTwoFreeformZAView, double>(ref hardOut, value); }
        }

        #endregion

        static LensTwoFreeformZAView()
        {
            InitRules();
        }

        internal LensTwoFreeformZAView(LensTwoFreeformZAParams ZA, MaterialView material)
        {
            modelObject = ZA;
            rInn = ZA.RInn;
            rOut = ZA.ROut;

            skirtHeight = ZA.SkirtHeight;
            gamma0 = ZA.Gamma0;
            hardOut = ZA.HardOut;

            Wavelength = ZA.Wavelength;
            Material = material;
            SubscribeToModelObject<LensTwoFreeformZAParams, LensTwoFreeformZAView>(ZA);
        }

        public static LensTwoFreeformZAView DefaultZAView()
        {
            var material = ApplicationPreferences.MaterialsBase.DefaultMaterial();
            var defaultZA =
                new LensTwoFreeformZAParams(material.modelObject);
            return new LensTwoFreeformZAView(defaultZA, material);
        }
        //public
        public override OpticalSceneObject ComputeOpticalElement(OpticalEntityView source,
                                                         RequiredSimulationResultView reqSim)
        {
            var inten = ((SolidSourceView)source.OEOrigin).IntensityModel.IntensityDistribution;
            var cs = source.CS;
            var lens = new LensTwoFreeform(modelObject, (IntensityAxisym)inten.ModelObject,
                                                 reqSim.ModelObject, cs.ModelObject);
            return new OpticalSceneObject(lens);
        }
        //private
        private static void InitRules()
        {
            #region validation
            RuleBook.ValidationRules.Add(new DelegateValidationRule<LensTwoFreeformZAView>(
                "RInn",
                "RInnOutOfRange",
                x => (Epsilon.CompareNumerics(x.RInn, 0) > 0 && Epsilon.CompareNumerics(x.RInn, x.ROut) < 0))
                );

            RuleBook.ValidationRules.Add(new DelegateValidationRule<LensTwoFreeformZAView>(
                "ROut",
                "ROutOutOfRange",
                x => (Epsilon.CompareNumerics(x.ROut, 0) > 0 && Epsilon.CompareNumerics(x.RInn, x.ROut) < 0)));

            RuleBook.ValidationRules.Add(new DelegateValidationRule<LensTwoFreeformZAView>(
                "SkirtHeight",
                "SkirtHeightLessThanZero",
                x => (Epsilon.CompareNumerics(x.SkirtHeight, 0) >= 0)));

            RuleBook.ValidationRules.Add(new DelegateValidationRule<LensTwoFreeformZAView>(
                "Gamma0",
                "Gamma0LessThanZero",
                x => (Epsilon.CompareNumerics(x.Gamma0, 0) >= 0)));

            RuleBook.ValidationRules.Add(new DelegateValidationRule<LensTwoFreeformZAView>(
                "HardOut",
                "HardOutOfRange",
                x => ((Epsilon.CompareNumerics(x.HardOut, 0) >= 0) && Epsilon.CompareNumerics(x.HardOut, 1) <= 0)
                ));

            #endregion

            #region bind
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

            RuleBook.BindRules.Add(new BindRule(
                "HardOut",
                "HardOut",
                x => x.modelObject.HardOut = x.HardOut,
                x => x.hardOut = x.modelObject.HardOut
                ));

            #endregion
        }

    }
}
