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
    using RuleBook = RuleBook<LensAxisymTIRAsphericZAView>;
    using ValidationRule = DelegateValidationRule<LensAxisymTIRAsphericZAView>;
    using BindRule = PropertyBindRule<LensAxisymTIRAsphericZAView>;

    public class LensAxisymTIRAsphericZAView : SolidLensZAParamsView
    {
        private LensAxisymTIRAsphericZAParams modelObject;

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
                SetProperty <LensAxisymTIRAsphericZAView,double>(ref r0, value);
            }
        }
        public double PsiBound
        {
            get { return UnitsConverter.RadToDeg(psiBound); }
            set
            {
                SetProperty<LensAxisymTIRAsphericZAView, double>
                    (ref psiBound, UnitsConverter.DegToRad(value));
            }
        }
        public double InnlatConicAng
        {
            get { return UnitsConverter.RadToDeg(innlatConicAng); }
            set
            {
                SetProperty<LensAxisymTIRAsphericZAView, double>(ref innlatConicAng, UnitsConverter.DegToRad(value));
            }
        }
        public double NeckHeight
        {
            get { return neckHeight; }
            set
            {
                SetProperty<LensAxisymTIRAsphericZAView, double>(ref neckHeight, value);
            }
        }
        public double SkirtWidth
        {
            get { return skirtWidth; }
            set
            {
                SetProperty<LensAxisymTIRAsphericZAView, double>(ref skirtWidth, value);
            }
        }
        public bool IsConvex
        {
            get { return isConvex; }
            set
            {
                SetProperty<LensAxisymTIRAsphericZAView, bool>(ref isConvex, value);
            }
        }
        # endregion

        //Constructors
        static LensAxisymTIRAsphericZAView()
        {
            InitRules();
        }

        internal LensAxisymTIRAsphericZAView(LensAxisymTIRAsphericZAParams ZA, MaterialView material)
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
            SubscribeToModelObject<LensAxisymTIRAsphericZAParams, LensAxisymTIRAsphericZAView>(ZA);
        }
        public static LensAxisymTIRAsphericZAView DefaultZAView()
        {
            var material = ApplicationPreferences.MaterialsBase.DefaultMaterial();
            var defaultZA =
                new LensAxisymTIRAsphericZAParams(material.modelObject);
            return new LensAxisymTIRAsphericZAView(defaultZA, material);
        }
        //public
        public override OpticalSceneObject ComputeOpticalElement(OpticalEntityView source,
                                                         RequiredSimulationResultView reqSim)
        {
            var inten = ((SolidSourceView)source.OEOrigin).IntensityModel.IntensityDistribution;
            var cs = source.CS;
            var lens = new LensAxisymTIRAspheric(modelObject, (IntensityAxisym)inten.ModelObject,
                                                 reqSim.ModelObject, cs.ModelObject);
            return new OpticalSceneObject(lens);
        }

        //Private methods
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
                    Epsilon.CompareNumerics(x.PsiBound, Math.PI / 3) < 0)));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "InnlatConicAng",
                "InnlatConicAngOutOfRange",
                x => (Epsilon.CompareNumerics(x.InnlatConicAng, 0) > 0 &&
                    Epsilon.CompareNumerics(x.InnlatConicAng, Math.PI / 4) < 0)));
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
