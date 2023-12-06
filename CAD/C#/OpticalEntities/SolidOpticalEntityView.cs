using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;
using ViewModelCore.MaterialProperties;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{
    public class SolidOpticalEntityView : OpticalEntityView, MaterialKeeper
    {
        // Fields
        private SolidOpticalEntity modelObject
        {
            get { return (SolidOpticalEntity)base.ModelObject; }
        }
        private MaterialView material;

        // Properties
        public MaterialView Material
        {
            get
            {
                if(material == null)
                    material = MaterialView.CreateByMaterial(modelObject.Material);
                return material;
            }
            set { SetProperty(ref material, value); }
        }

        public SolidOpticalEntityView(SolidOpticalEntity solid)
            :base(solid)
        {
            SubscribeToModelObject<SolidOpticalEntity, SolidOpticalEntityView>(solid);
        }

        internal new static SolidOpticalEntityView Wrap(OpticalEntity modelObj)
        {
            if (modelObj is LensAxisymFresnelTIR)
                return new LensAxisymFresnelTIRView((LensAxisymFresnelTIR)modelObj);
            if (modelObj is LensAxisymTIRAspheric)
                return new LensAxisymTIRAsphericView((LensAxisymTIRAspheric)modelObj);
            if (modelObj is LensAxisymTIRFreeform)
                return new LensAxisymTIRFreeformView((LensAxisymTIRFreeform)modelObj);
            if (modelObj is LensAxisymTIRUpperFlat)
                return new LensAxisymTIRUpperFlatView((LensAxisymTIRUpperFlat)modelObj);
            if (modelObj is LensAxisymTwoAspheric)
                return new LensAxisymTwoAsphericView((LensAxisymTwoAspheric)modelObj);
            if (modelObj is LensFreeformNoDome)
                return new LensFreeformNoDomeView((LensFreeformNoDome)modelObj);
            if (modelObj is LensFreeformSphDome)
                return new LensFreeformSphDomeView((LensFreeformSphDome)modelObj);
            if (modelObj is LensTwoFreeform)
                return new LensTwoFreeformView((LensTwoFreeform)modelObj);
            if (modelObj is SolidOpticalEntity)
                return new SolidOpticalEntityView((SolidOpticalEntity)modelObj);

            return null;
        }

        static SolidOpticalEntityView()
        {
            InitRules();
        }

        private static void InitRules()
        {
            RuleBook < SolidOpticalEntityView>.BindRules.Add(new PropertyBindRule<SolidOpticalEntityView>
            (
                "Material",
                "Material",
                x => x.modelObject.Material = x.material.modelObject,
                x => x.material = null
            ));
        }
    }
}
