using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;
using ViewModelCore.MaterialProperties;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{
    public class MirrorFreeformView : OpticalEntityView
    {
         // Field
        private OpticalSurfacePropertiesView surfProperties;
        public static int numSurfaces = 1;

        // Properties
        public OpticalSurfacePropertiesView SurfProperties
        {
            get
            {
                if (surfProperties == null)
                    surfProperties = new OpticalSurfacePropertiesView(((MirrorFreeform)ModelObject).SurfProperties);
                return surfProperties;
            }
        }
        public int NumSurfaces
        {
            get { return numSurfaces; }
        }

        public MirrorFreeformView(MirrorFreeform mirror)
            : base(mirror)
        {
            SubscribeToModelObject<MirrorFreeform, MirrorFreeformView>(mirror);
        }

        static MirrorFreeformView()
        {
            InitRules();
        }

        private static void InitRules()
        {
            RuleBook<MirrorFreeformView>.BindRules.Add(new PropertyBindRule<MirrorFreeformView>
                (
                "SurfProperties",
                "SurfProperties",
                null,
                x => x.surfProperties = null
            ));
        }
    }
}
