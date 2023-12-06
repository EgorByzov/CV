using System;
using MathKernel;
using ObjectModel.Common;
using ObjectModel.LightDistributions;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;
using ViewModelCore.MaterialProperties;
using ViewModelCore.OpticalSceneEntities.OriginEntities;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{
    using RuleBook = RuleBook<MirrorFreeformZAView>;
    using ValidationRule = DelegateValidationRule<MirrorFreeformZAView>;
    using BindRule = PropertyBindRule<MirrorFreeformZAView>;

    public class MirrorFreeformZAView : ZAParamsView
    {
        private double r0, alphaMin;
        private OpticalSurfacePropertiesView surfProperties;
        public MirrorFreeformZAParams modelObject;

        public double R0
        {
            get { return r0; }
            set
            {
                SetProperty<MirrorFreeformZAView, double>(ref r0, value);
            }
        }

        public double AlphaMin
        {
            get { return UnitsConverter.RadToDeg(alphaMin); }
            set
            {
                SetProperty<MirrorFreeformZAView, double>(ref alphaMin, UnitsConverter.DegToRad(value));
            }
        }

        public OpticalSurfacePropertiesView SurfProperties
        {
            get
            {
                if(surfProperties == null)
                    surfProperties = new OpticalSurfacePropertiesView(modelObject.SurfProperties);
                return surfProperties;
            }
            set
            {
                SetProperty<MirrorFreeformZAView, OpticalSurfacePropertiesView>(ref surfProperties, value);
            }
        }
        public MirrorFreeformZAView(MirrorFreeformZAParams ZA)
        {
            modelObject = ZA;
            r0 = ZA.R0;
            alphaMin = ZA.AlphaMin;
            SubscribeToModelObject<MirrorFreeformZAParams, MirrorFreeformZAView>(ZA);
        }
        public static MirrorFreeformZAView DefaultZAView()
        {
            return new MirrorFreeformZAView(new MirrorFreeformZAParams());
        }

        static MirrorFreeformZAView()
        {
            InitRules();
        }

        //Private methods

        private static void InitRules()
        {
            RuleBook.ValidationRules.Add(new ValidationRule(
                "R0", 
                "R0OutOfRange",
                x => Epsilon.CompareNumerics(x.R0, 0) > 0));

            RuleBook.ValidationRules.Add(new ValidationRule(
                "AlphaMin", 
                "AlphaMinOutOfRange",
                x => ((Epsilon.CompareNumerics(x.alphaMin, Math.PI / 36) > 0)
                    && (Epsilon.CompareNumerics(x.alphaMin, Math.PI / 2.25) < 0))));

            RuleBook.BindRules.Add(new BindRule(
                "R0",
                "R0",
                x => x.modelObject.R0 = x.R0,
                x => x.r0 = x.modelObject.R0
                ));
            RuleBook.BindRules.Add(new BindRule(
                "AlphaMin",
                "AlphaMin",
                x => x.modelObject.AlphaMin = Math.PI/2 - x.alphaMin,
                x => x.alphaMin = Math.PI / 2 - x.modelObject.AlphaMin
                ));
            RuleBook.BindRules.Add(new BindRule(
              "SurfProperties",
              "SurfProperties",
              x => x.modelObject.SurfProperties = x.SurfProperties.modelObject,
              x => x.surfProperties = null
              ));
        }


        public override OpticalSceneObject ComputeOpticalElement(OpticalEntityView source, Light_Distributions.RequiredSimulationResultView reqSim)
        {
            var inten = ((SolidSourceView)source.OEOrigin).IntensityModel.IntensityDistribution;
            var cs = source.CS;
            var lens = MirrorFreeform.CreateMirrorFreeformSurface(modelObject, (IntensityAxisym)inten.ModelObject,
                                                 reqSim.ModelObject, cs.ModelObject);
            return new OpticalSceneObject(lens);
        }
    }
}
