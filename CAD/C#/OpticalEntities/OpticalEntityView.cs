using System;
using System.Collections.Generic;
using System.Linq;
using Common.Memento;
using Newtonsoft.Json;
using ObjectModel.AnalyticalEngine;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ObjectModel.OpticalSceneEntities.OriginEntities;
using ViewModelCore.Common;
using ViewModelCore.OpticalSceneEntities.OriginEntities;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{

    public class OpticalEntityView : NotifyDataErrorInfo
    {
        private static Dictionary<Type, ZAParamsView> LensTypeModel2DefaultZA = new Dictionary<Type, ZAParamsView>
        {
            {typeof (LensAxisymTIRFreeform), LensAxisymTIRFreeformZAView.DefaultZAView()},
            {typeof (LensAxisymFresnelTIR), LensAxisymFresnelTIRZAView.DefaultZAView()},
            {typeof (LensAxisymTIRAspheric), LensAxisymTIRAsphericZAView.DefaultZAView()},
            {typeof (LensAxisymTIRUpperFlat),LensAxisymTIRUpperFlatZAView.DefaultZAView()},
            {typeof (LensAxisymTwoAspheric), LensAxisymTwoAsphericZAView.DefaultZAView()},
            {typeof (LensFreeformNoDome),LensFreeformNoDomeZAView.DefaultZAView()},
            {typeof (LensFreeformSphDome), LensFreeformSphDomeZAView.DefaultZAView()},
            {typeof (LensTwoFreeform), LensTwoFreeformZAView.DefaultZAView()},
            {typeof (MirrorFreeform), MirrorFreeformZAView.DefaultZAView()}
        };

        public static ZAParamsView GetDefaultZAByPrototype(string prototype)
        {

            var modelType = AdvisorSelectionOpticalElement.LensPrototype.FirstOrDefault(x => x.Value == prototype).Key;
            return LensTypeModel2DefaultZA[modelType];
        }
        #region fields
        [field: NonSerialized] CoordinateSystemView cs;
        [field: NonSerialized][field:ReferenceMemorizableAttribute] OpticalEntityView parent;
        [field: NonSerialized] List<OpticalEntityView> children;
        [field: NonSerialized] protected OpticalEntityOriginView oeOrigin;
        bool isRaytraceOn;
        bool isOptimizationOn;
        private OpticalEntity modelObject;
        #endregion

        #region Constructors

        static OpticalEntityView()
        {
            InitRules();
        }
        [JsonConstructor]
        protected OpticalEntityView(OpticalEntity modelObject)
        {
            this.modelObject = modelObject;

            isRaytraceOn = modelObject.IsRaytraceOn;
            isOptimizationOn = modelObject.IsOptimizationOn;
            SubscribeToModelObject<OpticalEntity, OpticalEntityView>(modelObject);
        }

        internal static OpticalEntityView Wrap(OpticalEntity modelObj)
        {
            if (modelObj is Registrator)
                return new RegistratorView((Registrator)modelObj);
            if (modelObj is MirrorFreeform)
                return new MirrorFreeformView((MirrorFreeform)modelObj);
            if (modelObj is SolidOpticalEntity)
                return SolidOpticalEntityView.Wrap(modelObj);

            return new OpticalEntityView(modelObj);
        }

        #endregion

        #region Properties

        protected OpticalEntity ModelObject
        {
            get { return modelObject; }
        }

        public OpticalEntityOriginView OEOrigin
        {
            get
            {
                if (oeOrigin == null) oeOrigin = OpticalEntityOriginView.Wrap(ModelObject.OEOrigin);
                return oeOrigin;
            }
        }
        public CoordinateSystemView CS
        {
            get
            {
                if (cs == null) cs = new CoordinateSystemView(ModelObject.CS, ModelObject.ParentCSInGlobal) { IsEditEnabled = IsCSEditEnabled };
                return cs;
            }
        }
        public OpticalEntityView Parent
        {
            get
            {
                if(parent == null) parent = new OpticalEntityView(ModelObject.Parent);
                return parent;
            }
        }
        public List<OpticalEntityView> Children
        {
            get
            {
                if (children != null) return children;
                children = new List<OpticalEntityView>();
                foreach (var child in ModelObject.Children)
                {
                    children.Add(new OpticalEntityView(child));
                }
                return children;
            }
        }
        public bool IsRaytraceOn
        {
            get { return isRaytraceOn; }
            set { SetProperty(ref isRaytraceOn, value); }
        }
        public bool IsOptimizationOn
        {
            get { return isOptimizationOn; }
            set { SetProperty(ref isOptimizationOn, value); }
        }
        public bool IsCSEditEnabled
        {
            get
            {
                if (ModelObject == null || ModelObject.Parent == null) return true;
                return !ModelObject.Parent.IsSolidPart;
            }
        }
        #endregion

        #region Private
        private static void InitRules()
        {
            RuleBook<OpticalEntityView>.BindRules.Add(new PropertyBindRule<OpticalEntityView>
              (
              "IsRaytraceOn",
              "IsRaytraceOn",
              x => x.ModelObject.IsRaytraceOn = x.IsRaytraceOn,
              x => x.IsRaytraceOn = x.ModelObject.IsRaytraceOn
              ));
            RuleBook<OpticalEntityView>.BindRules.Add(new PropertyBindRule<OpticalEntityView>
              (
              "IsOptimizationOn",
              "IsOptimizationOn",
              x => x.ModelObject.IsOptimizationOn = x.IsOptimizationOn,
              x => x.isOptimizationOn = x.ModelObject.IsOptimizationOn
              ));
            RuleBook<OpticalEntityView>.BindRules.Add(new PropertyBindRule<OpticalEntityView>
              (
              "CS",
              "CS",
              null,
              x => x.cs = null
              ));
            RuleBook<OpticalEntityView>.BindRules.Add(new PropertyBindRule<OpticalEntityView>
              (
              "Parent",
              "Parent",
              null,
              x => x.parent = null
              ));
            RuleBook<OpticalEntityView>.BindRules.Add(new PropertyBindRule<OpticalEntityView>
              (
              "Children",
              "Children",
              null,
              x => x.children = null
              ));
            RuleBook<OpticalEntityView>.BindRules.Add(new PropertyBindRule<OpticalEntityView>
              (
              "IsSolidPart",
              "IsSolidPart",
              null,
              x => x.OnPropertyChanged("IsSolidPart")
              ));
        }

        #endregion
    }
}
