using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathKernel;
using ObjectModel.OpticalSceneEntities;
using ViewModelCore.Common;

namespace ViewModelCore.OpticalSceneEntities
{
    using RuleBook = RuleBook<CylinderView>;
    using ValidationRule = DelegateValidationRule<CylinderView>;
    using BindRule = PropertyBindRule<CylinderView>;

    public class CylinderView : EmittableSolidView
    {
        // Fields
        private Cylinder modelObject;
        private double radius;
        private double height;


        // Properties
        public override EmittableSolidType Type
        {
            get { return EmittableSolidType.Cylinder; }
        }

        internal override EmittableSolid ModelObject
        {
            get { return modelObject; }
        }
        public double Radius
        {
            get { return radius; }
            set { SetProperty<CylinderView, double>(ref radius, value); }
        }
        public double Height
        {
            get { return height; }
            set { SetProperty<CylinderView, double>(ref height, value); }
        }

        internal CylinderView(Cylinder cylinder)
        {
            modelObject = cylinder;
            radius = cylinder.Radius;
            height = cylinder.Height;

            SubscribeToModelObject<Cylinder , CylinderView>(cylinder);
        }
        static CylinderView()
        {
            InitRules();
        }

        private static void InitRules()
        {
            RuleBook.ValidationRules.Add(new ValidationRule(
              "Radius",
              "RadiusOutOfRange",
              x => Epsilon.CompareNumerics(x.Radius, 0) >= 0)
              );
            RuleBook.ValidationRules.Add(new ValidationRule(
                "Height",
                "HeightOutOfRange",
                x => Epsilon.CompareNumerics(x.Height, 0) >= 0)
                );

            RuleBook.BindRules.Add(new BindRule
               (
               "Radius",
               "Radius",
               x => x.modelObject.Radius = x.Radius,
               x => x.radius = x.modelObject.Radius
               ));
            RuleBook.BindRules.Add(new BindRule
               (
               "Height",
               "Height",
               x => x.modelObject.Height = x.Height,
               x => x.height = x.modelObject.Height
               ));
        }

    }
}
