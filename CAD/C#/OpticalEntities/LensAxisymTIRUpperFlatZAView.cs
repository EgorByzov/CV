using System;
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
    using RuleBook = RuleBook<LensAxisymTIRUpperFlatZAView>;
    using ValidationRule = DelegateValidationRule<LensAxisymTIRUpperFlatZAView>;
    using BindRule = PropertyBindRule<LensAxisymTIRUpperFlatZAView>;

    public class LensAxisymTIRUpperFlatZAView : SolidLensZAParamsView
    {
        private LensAxisymTIRUpperFlatZAParams modelObject;
        public override SolidLensZAParams ModelObject
        {
            get { return modelObject; }
        }
        private double r0, psiBound, innlatConicAng,
                       neckHeight, skirtWidth;


        # region Properties
        public double R0
        {
            get { return r0; }
            set
            {
                SetProperty < LensAxisymTIRUpperFlatZAView, double>(ref r0, value);
            }
        }
        public double PsiBound
        {
            get { return psiBound; }
            set
            {
                SetProperty<LensAxisymTIRUpperFlatZAView, double>
                    (ref psiBound, value);
            }
        }
        public double InnlatConicAng
        {
            get { return innlatConicAng; }
            set
            {
                SetProperty<LensAxisymTIRUpperFlatZAView, double>
                    (ref innlatConicAng, value);
            }
        }
        public double NeckHeight
        {
            get { return neckHeight; }
            set
            {
                SetProperty<LensAxisymTIRUpperFlatZAView, double>(ref neckHeight, value);
            }
        }
        public double SkirtWidth
        {
            get { return skirtWidth; }
            set
            {
                SetProperty<LensAxisymTIRUpperFlatZAView, double>(ref skirtWidth, value);
            }
        }
        # endregion

        // Constructors
        static LensAxisymTIRUpperFlatZAView()
        {
            InitRules();
        }

        internal LensAxisymTIRUpperFlatZAView(LensAxisymTIRUpperFlatZAParams ZA, MaterialView material)
        {
            modelObject = ZA;
            r0 = ZA.R0;
            psiBound = ZA.PsiBound;
            innlatConicAng = ZA.InnlatConicAng;
            neckHeight = ZA.NeckHeight;
            skirtWidth = ZA.SkirtWidth;
            Wavelength = ZA.Wavelength;
            Material = material;

            SubscribeToModelObject<LensAxisymTIRUpperFlatZAParams, LensAxisymTIRUpperFlatZAView>(ZA);
        }


        public static LensAxisymTIRUpperFlatZAView DefaultZAView()
        {
            var material = ApplicationPreferences.MaterialsBase.DefaultMaterial();
            var defaultZA =
                new LensAxisymTIRUpperFlatZAParams(material.modelObject);
            return new LensAxisymTIRUpperFlatZAView(defaultZA, material);
        }

        //public
        public override OpticalSceneObject ComputeOpticalElement(OpticalEntityView source,
                                                         RequiredSimulationResultView reqSim)
        {
            var inten = ((SolidSourceView)source.OEOrigin).IntensityModel.IntensityDistribution;
            var cs = source.CS;
            var lens = new LensAxisymTIRUpperFlat(modelObject, (IntensityAxisym)inten.ModelObject,
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
                    Epsilon.CompareNumerics(x.PsiBound, Math.PI / 3) <= 0)));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "InnlatConicAng",
                "InnlatConicAngOutOfRange",
                x => (Epsilon.CompareNumerics(x.PsiBound, 0) > 0 &&
                    Epsilon.CompareNumerics(x.PsiBound, Math.PI / 4) <= 0)));
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
                x => x.modelObject.PsiBound = x.PsiBound,
                x => x.PsiBound = x.modelObject.PsiBound
                ));
            RuleBook.BindRules.Add(new BindRule(
                "InnlatConicAng",
                "InnlatConicAng",
                x => x.modelObject.InnlatConicAng = x.InnlatConicAng,
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
