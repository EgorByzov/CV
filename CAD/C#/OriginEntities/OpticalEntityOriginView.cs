using System;
using ObjectModel.OpticalSceneEntities.OriginEntities;
using ObjectModel.OpticalSceneEntities.OriginEntities.Surfaces;
using ViewModelCore.Common;
using ViewModelCore.OpticalSceneEntities.OriginEntities.Surfaces;

namespace ViewModelCore.OpticalSceneEntities.OriginEntities
{
    public class OpticalEntityOriginView : NotifyDataErrorInfo
    {
        // Fields
        public OpticalEntityOrigin ModelObject { get; protected set; }

        protected OpticalEntityOriginView(OpticalEntityOrigin oe)
        {
            ModelObject = oe;
            if (ModelObject != null)
            {
                SubscribeToModelObject<OpticalEntityOrigin, OpticalEntityOriginView>(oe);
            }
        }

        static OpticalEntityOriginView()
        {
            InitRules();
        }

        private static void InitRules()
        {
        }

        internal static OpticalEntityOriginView Wrap(OpticalEntityOrigin obj)
        {
            if (obj is SolidSource)
                return new SolidSourceView((SolidSource)obj);

            if (obj is Surface)
                return SurfaceView.Wrap((Surface) obj);

            return null;
        }
    }
}
