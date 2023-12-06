using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Markup;
using MathKernel;
using ObjectModel.OpticalSceneEntities;
using ViewModelCore.Common;

namespace ViewModelCore.OpticalSceneEntities
{
    using RuleBook = RuleBook<ParallelepipedView>;
    using ValidationRule = DelegateValidationRule<ParallelepipedView>;
    using BindRule = PropertyBindRule<ParallelepipedView>;

    public class ParallelepipedView : EmittableSolidView
    {
        // Fields
        private Parallelepiped modelObject;
        private double xSize;
        private double ySize;
        private double zSize;

        // Properties
        public override EmittableSolidType Type
        {
            get { return EmittableSolidType.Parallelepiped; }
        }

        internal override EmittableSolid ModelObject
        {
            get { return modelObject; }
        }
        public double XSize
        {
            get { return xSize; }
            set { SetProperty<ParallelepipedView, double>(ref xSize, value); }
        }
        public double YSize
        {
            get { return ySize; }
            set { SetProperty<ParallelepipedView, double>(ref ySize, value); }
        }
        public double ZSize
        {
            get { return zSize; }
            set { SetProperty<ParallelepipedView, double>(ref zSize, value); }
        }

        internal ParallelepipedView(Parallelepiped parallelepiped)
        {
            modelObject = parallelepiped;
            xSize = parallelepiped.XSize;
            ySize = parallelepiped.YSize;
           zSize = parallelepiped.ZSize;

            SubscribeToModelObject<Parallelepiped, ParallelepipedView>(parallelepiped);
        }

        static ParallelepipedView()
        {
            InitRules();
        }
        public static void InitRules()
        {
            RuleBook.ValidationRules.Add(new ValidationRule(
                "XSize",
                "XSizeOutOfRange",
                x => Epsilon.CompareNumerics(x.XSize, 0) >= 0)
                );
            RuleBook.ValidationRules.Add(new ValidationRule(
                "YSize",
                "YSizeOutOfRange",
                x => Epsilon.CompareNumerics(x.YSize, 0) >= 0)
                );
            RuleBook.ValidationRules.Add(new ValidationRule(
                "ZSize",
                "ZSizeOutOfRange",
                x => Epsilon.CompareNumerics(x.ZSize, 0) >= 0)
                );

            RuleBook.BindRules.Add(new BindRule
                (
                "XSize",
                "XSize",
                x => x.modelObject.XSize = x.XSize,
                x => x.xSize = x.modelObject.XSize
                ));
            RuleBook.BindRules.Add(new BindRule
                (
                "YSize",
                "YSize",
                x => x.modelObject.YSize = x.YSize,
                x => x.ySize = x.modelObject.YSize
                ));
            RuleBook.BindRules.Add(new BindRule
                (
                "ZSize",
                "ZSize",
                x => x.modelObject.ZSize = x.ZSize,
                x => x.zSize = x.modelObject.ZSize
                ));
        }

        
    }
}
