using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ObjectModel.OpticalSceneEntities.OriginEntities.Surfaces;
using ViewModelCore.Common;
using ViewModelCore.MaterialProperties;

namespace ViewModelCore.OpticalSceneEntities.OriginEntities.Surfaces
{
    public class SurfaceView : OpticalEntityOriginView
    {
        // fields
        private Surface modelObject
        {
            get { return (Surface) base.ModelObject; }
        }
        private OpticalSurfacePropertiesView opticalProperties;

        // Properties
        public OpticalSurfacePropertiesView OpticalProperties
        {
            get
            {
                if(opticalProperties == null) opticalProperties = new OpticalSurfacePropertiesView(modelObject.OpticalProperties);
                return opticalProperties;
            }
        }

        internal SurfaceView(Surface surface)
            :base(surface)
        {

            SubscribeToModelObject<Surface, SurfaceView>(surface);
        }

        static SurfaceView()
        {
            InitRules();
        }

        private static void InitRules()
        {
            RuleBook<SurfaceView>.BindRules.Add(new PropertyBindRule<SurfaceView>
                (
                "OpticalProperties",
                "OpticalProperties",
                null,
                x => x.opticalProperties = null
                ));
        }

        internal static SurfaceView Wrap(Surface obj)
        {
            if (obj is AxisymBezierSurface)
                return new AxisymBezierSurfaceView((AxisymBezierSurface)obj);
            if (obj is AxisymSphereSegment)
                return new AxisymSphereSegmentView((AxisymSphereSegment)obj);
            if (obj is FreeformSurface)
                return new FreeformSurfaceView((FreeformSurface)obj);
            if (obj is RectangularSurface)
                return new RectangularSurfaceView((RectangularSurface)obj);

            return new SurfaceView((Surface)obj);
        }
    }
}
