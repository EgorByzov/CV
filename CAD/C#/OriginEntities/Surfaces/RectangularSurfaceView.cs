using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathKernel;
using ObjectModel.OpticalSceneEntities.OriginEntities.Surfaces;
using ViewModelCore.Common;

namespace ViewModelCore.OpticalSceneEntities.OriginEntities.Surfaces
{
    using RuleBook = RuleBook<RectangularSurfaceView>;
    using ValidationRule = DelegateValidationRule<RectangularSurfaceView>;
    using BindRule = PropertyBindRule<RectangularSurfaceView>;

    public class RectangularSurfaceView : SurfaceView
    {
        private RectangularSurface modelObject
        {
            get { return (RectangularSurface) base.ModelObject; }
        }

        private double xSize, ySize;

        // Properties
        public double XSize
        {
            get { return xSize; }
            set { SetProperty<RectangularSurfaceView, double>(ref xSize, value); }
        }
        public double YSize
        {
            get { return ySize; }
            set { SetProperty<RectangularSurfaceView, double>(ref ySize, value); }
        }

        public RectangularSurfaceView(RectangularSurface recSurface)
            :base(recSurface)
        {
            // Init fields
            xSize = recSurface.XSize;
            ySize = recSurface.YSize;

            SubscribeToModelObject<RectangularSurface, RectangularSurfaceView>(recSurface);
        }

        private static void InitRules()
        {
            RuleBook.ValidationRules.Add(new ValidationRule(
                "XSize",
                "XSizeOutOfRange",
                x => Epsilon.CompareNumerics(x.XSize, 0) > 0)
                );
            RuleBook.BindRules.Add(new BindRule(
                "XSize",
                "XSize",
                x => x.modelObject.XSize = x.XSize,
                x => x.xSize = x.modelObject.XSize
            ));
            RuleBook.ValidationRules.Add(new ValidationRule(
                "YSize",
                "YSizeOutOfRange",
                x => Epsilon.CompareNumerics(x.YSize, 0) > 0)
                );
            RuleBook.BindRules.Add(new BindRule(
                "YSize",
                "YSize",
                x => x.modelObject.YSize = x.YSize,
                x => x.ySize = x.modelObject.YSize
            ));

        }

        static RectangularSurfaceView()
        {
            InitRules();
        }

    }
}
