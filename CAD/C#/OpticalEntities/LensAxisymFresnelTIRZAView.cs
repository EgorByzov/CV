using MathKernel;
using MathKernel.VectorAlgebra;
using ObjectModel.Common;
using ObjectModel.LightDistributions;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;
using ViewModelCore.Light_Distributions;
using ViewModelCore.MaterialProperties;
using ViewModelCore.OpticalSceneEntities.OriginEntities;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{
    using RuleBook = RuleBook<LensAxisymFresnelTIRZAView>;
    using ValidationRule = DelegateValidationRule<LensAxisymFresnelTIRZAView>;
    using BindRule = PropertyBindRule<LensAxisymFresnelTIRZAView>;

    public class LensAxisymFresnelTIRZAView : SolidLensZAParamsView
    {
        private LensAxisymFresnelTIRZAParams modelObject;

        public override SolidLensZAParams ModelObject
        {
            get { return modelObject; }
        }

        private double h0, hRelief, hBase, alphaInc,
                       rRefract, rTIR, rMax;

        # region Properties
        public double H0
        {
            get { return h0; }
            set
            {
                SetProperty<LensAxisymTwoAsphericZAView,double>(ref h0, value);
            }
        }
        public double HRelief
        {
            get { return hRelief; }
            set
            {
                SetProperty<LensAxisymTwoAsphericZAView, double>(ref hRelief, value);
            }
        }
        public double HBase
        {
            get { return hBase; }
            set
            {
                SetProperty<LensAxisymTwoAsphericZAView, double>(ref hBase, value);
            }
        }
        public double AlphaInc
        {
            get { return UnitsConverter.RadToDeg(alphaInc); }
            set
            {
                SetProperty<LensAxisymTwoAsphericZAView, double>(ref alphaInc, UnitsConverter.DegToRad(value));
            }
        }
        public double RRefract
        {
            get { return rRefract; }
            set
            {
                SetProperty<LensAxisymTwoAsphericZAView, double>(ref rRefract, value);
            }
        }
        public double RTIR
        {
            get { return rTIR; }
            set
            {
                SetProperty<LensAxisymTwoAsphericZAView, double>(ref rTIR, value);
            }
        }
        public double RMax
        {
            get { return rMax; }
            set
            {
                SetProperty<LensAxisymTwoAsphericZAView, double>(ref rMax, value);
            }
        }
        # endregion

        //Constructors
        static LensAxisymFresnelTIRZAView()
        {
            InitRules();
        }

        internal LensAxisymFresnelTIRZAView(LensAxisymFresnelTIRZAParams ZA, MaterialView material)
        {
            modelObject = ZA;
            h0 = ZA.H0;
            hRelief = ZA.HRelief;
            hBase = ZA.HBase;
            alphaInc = ZA.AlphaInc;
            rRefract = ZA.RRefract;
            rTIR = ZA.RTIR;
            rMax = ZA.RMax;
            Wavelength = ZA.Wavelength;
            Material = material;
            SubscribeToModelObject<LensAxisymFresnelTIRZAParams, LensAxisymFresnelTIRZAView>(ZA);
        }
        public static LensAxisymFresnelTIRZAView DefaultZAView()
        {
            var material = ApplicationPreferences.MaterialsBase.DefaultMaterial();
            var defaultZA =
                new LensAxisymFresnelTIRZAParams(material.modelObject);
            return new LensAxisymFresnelTIRZAView(defaultZA, material);
        }
        //public
        public override OpticalSceneObject ComputeOpticalElement(OpticalEntityView source,
                                                         RequiredSimulationResultView reqSim)
        {
            var cs = source.CS;
            var lens = new LensAxisymFresnelTIR(modelObject, cs.ModelObject);
            return new OpticalSceneObject(lens);
        }

        //private
        private static void InitRules()
        {
            RuleBook.ValidationRules.Add(new ValidationRule(
                "H0",
                "H0OutOfRange",
                x => Epsilon.CompareNumerics(x.H0, 0) > 0));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "HRelief",
                "HReliefOutOfRange",
                x => Epsilon.CompareNumerics(x.HRelief, 0) > 0));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "HBase",
                "HBaseOutOfRange",
                x => Epsilon.CompareNumerics(x.HBase, 0) > 0));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "AlphaInc",
                "AlphaIncOutOfRange",
                x => Epsilon.CompareNumerics(x.AlphaInc, 0) > 0));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "RRefract",
                "RRefractOutOfRange",
                x => Epsilon.CompareNumerics(x.RRefract, 0) > 0));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "RTIR",
                "RTIROutOfRange",
                x => Epsilon.CompareNumerics(x.RTIR, 0) > 0));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "RMax",
                "RMaxOutOfRange",
                x => Epsilon.CompareNumerics(x.RMax, 0) > 0));

            //BindRule
            RuleBook.BindRules.Add(new BindRule(
               "H0",
               "H0",
                x => x.modelObject.H0 = x.H0,
                x => x.h0 = x.modelObject.H0));
            RuleBook.BindRules.Add(new BindRule(
                "HRelief",
                "HRelief",
                 x => x.modelObject.HRelief = x.HRelief,
                 x => x.hRelief = x.modelObject.HRelief));
            RuleBook.BindRules.Add(new BindRule(
                "HBase",
                "HBase",
                 x => x.modelObject.HBase = x.HBase,
                 x => x.hBase = x.modelObject.HBase));
            RuleBook.BindRules.Add(new BindRule(
                "AlphaInc",
                "AlphaInc",
                 x => x.modelObject.AlphaInc = x.AlphaInc,
                 x => x.alphaInc = x.modelObject.AlphaInc));
            RuleBook.BindRules.Add(new BindRule(
                "RRefract",
                "RRefract",
                x => x.modelObject.RRefract = x.RRefract,
                x => x.rRefract = x.modelObject.RRefract));
            RuleBook.BindRules.Add(new BindRule(
                "RTIR",
                "RTIR",
                x => x.modelObject.RTIR = x.RTIR,
                x => x.rTIR = x.modelObject.RTIR));
            RuleBook.BindRules.Add(new BindRule(
                "RMax",
                "RMax",
                x => x.modelObject.RMax = x.RMax,
                x => x.rMax = x.modelObject.RMax));
        }
    }
}
