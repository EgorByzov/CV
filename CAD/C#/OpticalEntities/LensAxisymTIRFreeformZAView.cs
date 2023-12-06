using MathKernel;
using ObjectModel.Common;
using ObjectModel.LightDistributions;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;
using ViewModelCore.Light_Distributions;
using ViewModelCore.MaterialProperties;
using ViewModelCore.OpticalSceneEntities.OriginEntities;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{
    using RuleBook = RuleBook<LensAxisymTIRFreeformZAView>;
    using ValidationRule = DelegateValidationRule<LensAxisymTIRFreeformZAView>;
    using BindRule = PropertyBindRule<LensAxisymTIRFreeformZAView>;
    
    public class LensAxisymTIRFreeformZAView : SolidLensZAParamsView
    {
        private LensAxisymTIRFreeformZAParams modelObject;
        public override SolidLensZAParams ModelObject
        {
            get { return modelObject; }
        }
        private double r0, psiBound, innlatConicAng,
                       neckHeight, skirtWidth;
        private bool isConvex;

        # region Properties
        public double R0
        {
            get { return r0; }
            set
            {
                SetProperty<LensAxisymTIRFreeformZAView, double>(ref r0, value);
            }
        }
        public double PsiBound
        {
            get { return UnitsConverter.RadToDeg(psiBound); }
            set
            {
                SetProperty<LensAxisymTIRFreeformZAView, double>
                    (ref psiBound, UnitsConverter.DegToRad(value));
            }
        }
        public double InnlatConicAng
        {
            get { return UnitsConverter.RadToDeg(innlatConicAng); }
            set
            {
                SetProperty<LensAxisymTIRFreeformZAView, double>
                    (ref innlatConicAng, UnitsConverter.DegToRad(value));
            }
        }
        public double NeckHeight
        {
            get { return neckHeight; }
            set
            {
                SetProperty<LensAxisymTIRFreeformZAView, double>(ref neckHeight, value);
            }
        }
        public double SkirtWidth
        {
            get { return skirtWidth; }
            set
            {
                SetProperty<LensAxisymTIRFreeformZAView, double>(ref skirtWidth, value);
            }
        }
        public bool IsConvex
        {
            get { return isConvex; }
            set
            {
                SetProperty<LensAxisymTIRFreeformZAView, bool>(ref isConvex, value);
            }
        }
        # endregion

        //Constructors
        static LensAxisymTIRFreeformZAView()
        {
            InitRules();
        }

        internal LensAxisymTIRFreeformZAView(LensAxisymTIRFreeformZAParams ZA,
            MaterialView material)
        {
            modelObject = ZA;
            r0 = ZA.R0;
            psiBound = ZA.PsiBound;
            innlatConicAng = ZA.InnlatConicAng;
            neckHeight = ZA.NeckHeight;
            skirtWidth = ZA.SkirtWidth;
            isConvex = ZA.IsConvex;
            Wavelength = ZA.Wavelength;
            Material = material;
            SubscribeToModelObject<LensAxisymTIRFreeformZAParams, LensAxisymTIRFreeformZAView>(ZA);
        }

        public static LensAxisymTIRFreeformZAView DefaultZAView()
        {
            var material = ApplicationPreferences.MaterialsBase.DefaultMaterial();
            var defaultZA =
                new LensAxisymTIRFreeformZAParams(material.modelObject);
            return new LensAxisymTIRFreeformZAView(defaultZA, material);
        }

        //public
        public override OpticalSceneObject ComputeOpticalElement(OpticalEntityView source,
                                                         RequiredSimulationResultView reqSim)
        {
            var sourceOrigin = (SolidSourceView) source.OEOrigin;
            var inten = sourceOrigin.IntensityModel.IntensityDistribution;
            var cs = source.CS;
            var lens = new LensAxisymTIRFreeform(modelObject, (IntensityAxisym)inten.ModelObject,
                                                 reqSim.ModelObject, cs.ModelObject);
            return new OpticalSceneObject(lens);
        }

        //private
        private static void InitRules()
        {
            RuleBook.ValidationRules.Add(new ValidationRule(
                "R0",
                "R0OutOfRange",
                x => Epsilon.CompareNumerics(x.R0, 0) > 0));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "PsiBound",
                "PsiBoundOutOfRange",
                x => (Epsilon.CompareNumerics(x.PsiBound, 0) > 0 &&
                    Epsilon.CompareNumerics(x.PsiBound, 60) < 0)));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "InnlatConicAng",
                "InnlatConicAngOutOfRange",
                x => (Epsilon.CompareNumerics(x.InnlatConicAng, 0) > 0 &&
                    Epsilon.CompareNumerics(x.InnlatConicAng, 45) < 0)));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "NeckHeight",
                "NeckHeightOutOfRange",
                x => Epsilon.CompareNumerics(x.NeckHeight, 0) > 0));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "SkirtWidth",
                "SkirtWidthOutOfRange",
                x => Epsilon.CompareNumerics(x.SkirtWidth, 0) > 0));

            RuleBook.BindRules.Add(new BindRule(
                "R0",
                "R0",
                x => x.modelObject.R0 = x.R0,
                x => x.r0 = x.modelObject.R0
                ));
            RuleBook.BindRules.Add(new BindRule(
                "PsiBound",
                "PsiBound",
                x => x.modelObject.PsiBound = x.psiBound,
                x => x.psiBound = x.modelObject.PsiBound
                ));
            RuleBook.BindRules.Add(new BindRule(
                "InnlatConicAng",
                "InnlatConicAng",
                x => x.modelObject.InnlatConicAng = x.innlatConicAng,
                x => x.innlatConicAng = x.modelObject.InnlatConicAng
                ));
            RuleBook.BindRules.Add(new BindRule(
                "NeckHeight",
                "NeckHeight",
                x => x.modelObject.NeckHeight = x.NeckHeight,
                x => x.neckHeight = x.modelObject.NeckHeight
                ));
            RuleBook.BindRules.Add(new BindRule(
                "SkirtWidth",
                "SkirtWidth",
                x => x.modelObject.SkirtWidth = x.SkirtWidth,
                x => x.skirtWidth = x.modelObject.SkirtWidth
                ));
        }


    }
}
