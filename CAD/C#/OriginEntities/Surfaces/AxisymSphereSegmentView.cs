using System;
using MathKernel;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ObjectModel.OpticalSceneEntities.OriginEntities.Surfaces;
using ViewModelCore.Common;
using ViewModelCore.OpticalSceneEntities.OpticalEntities;

namespace ViewModelCore.OpticalSceneEntities.OriginEntities.Surfaces
{
    using RuleBook = RuleBook<AxisymSphereSegmentView>;
    using ValidationRule = DelegateValidationRule<AxisymSphereSegmentView>;
    using BindRule = PropertyBindRule<AxisymSphereSegmentView>;

    public class AxisymSphereSegmentView :  SurfaceView
    {
        // fields
        private AxisymSphereSegment modelObject
        {
            get { return (AxisymSphereSegment) base.ModelObject; }
        }
        private double radius;
        private double psiMin;
        private double psiMax;

        // Properties
        public double Radius
        {
            get { return radius; }
            set { SetProperty <AxisymSphereSegmentView, double>(ref radius, value); }
        }
        public double PsiMin
        {
            get { return psiMin * 2 * 180 / Math.PI; }
            set { SetProperty<AxisymSphereSegmentView, double>(ref psiMin, value / 2 / 180 * Math.PI); }
        }
        public double PsiMax
        {
            get { return psiMax * 2 * 180 / Math.PI; }
            set { SetProperty<AxisymSphereSegmentView, double>(ref psiMax, value / 2 / 180 * Math.PI); }
        }

        public AxisymSphereSegmentView(AxisymSphereSegment sphere)
            : base(sphere)
        {
            radius = sphere.Radius;
            psiMin = sphere.PsiMin;
            psiMax = sphere.PsiMax;

            SubscribeToModelObject<AxisymSphereSegment, AxisymSphereSegmentView>(sphere);
        }

        static AxisymSphereSegmentView()
        {
            InitRules();
        }

        private static void InitRules()
        {
            RuleBook.ValidationRules.Add(new ValidationRule(
              "Radius",
              "RadiusOutOfRange",
              x => Epsilon.CompareNumerics(x.Radius, 0) > 0)
              );
            RuleBook.ValidationRules.Add(new ValidationRule(
                "PsiMax",
                "PsiMaxRange",
                x => ((Epsilon.CompareNumerics(x.psiMax, x.psiMin) >= 0) &&
                      (Epsilon.CompareNumerics(x.psiMax, 2 * Math.PI) <= 0)
                )));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "PsiMin",
                "PsiMinRange",
                x => ((Epsilon.CompareNumerics(x.psiMin, 0) >= 0) &&
                      (Epsilon.CompareNumerics(x.psiMin, x.psiMax) <= 0)
                )));
            RuleBook.BindRules.Add(new BindRule
               (
               "Radius",
               "Radius",
               x => x.modelObject.Radius = x.Radius,
               x => x.radius = x.modelObject.Radius
               ));
            RuleBook.BindRules.Add(new BindRule
               (
               "PsiMax",
               "PsiMax",
               x => x.modelObject.PsiMax = x.psiMax,
               x => x.psiMax = x.modelObject.PsiMax
               ));
            RuleBook.BindRules.Add(new BindRule
               (
               "PsiMin",
               "PsiMin",
               x => x.modelObject.PsiMin = x.psiMin,
               x => x.psiMin = x.modelObject.PsiMin
               ));

        }
    }
}
