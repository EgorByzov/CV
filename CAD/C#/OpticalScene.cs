using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Runtime.Serialization;
using Common;
using Newtonsoft.Json;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.OpticalSceneEntities.SceneView;

namespace ViewModelCore.OpticalSceneEntities
{
    [JsonObject(MemberSerialization.Fields)]
    public class OpticalScene
    {
        #region Fields

        private readonly OpticalEntity baseTree;
        private OpticalSceneObject sceneObjectTree;

        private readonly ObservableCollection<OpticalSceneObject> registrators;
        private readonly ObservableCollection<OpticalSceneObject> sources;

        [field:NonSerialized]
        private OpenGLOpticalScene scenePresenter;

        #endregion Fields

        [field: NonSerialized]
        public event EventHandler<ObjectChangedEventArgs> ObjectChanged;
        public enum EventIDs { Add, Replace, Delete }

        #region Properties

        public OpenGLOpticalScene ScenePresenter
        {
            get { return scenePresenter; }
        }

        internal OpticalEntity BaseTree
        {
            get { return baseTree; }
        }

        internal OpticalSceneObject SceneObjectTree
        {
            get
            {
                return sceneObjectTree ?? (sceneObjectTree = new OpticalSceneObject(baseTree));
            }
        }

        internal ObservableCollection<OpticalSceneObject> Registrators
        {
            get { return registrators; }
        }

        internal ObservableCollection<OpticalSceneObject> Sources
        {
            get { return sources; }
        }

        #endregion Properties

        internal OpticalScene()
        {
            baseTree = new OpticalEntity();
            registrators = new ObservableCollection<OpticalSceneObject>();
            sources = new ObservableCollection<OpticalSceneObject>();
            scenePresenter = new OpenGLOpticalScene(this);

            SceneObjectTree.RoutedChange += OnOSObjectChanged;
        }

        #region Methods

        public bool AddOpticalEntity(OpticalEntity newOE)
        {
            var newOSObject = new OpticalSceneObject(newOE);
            return AddOpticalSceneObject(newOSObject);
        }

        internal bool AddOpticalSceneObject(OpticalSceneObject newOSObject)
        {
            if (newOSObject.IsRegistrator) ScenePresenter.SetCurrentRoot(SceneObjectTree);

            var cRoot = ScenePresenter.CurrentRoot;
            while (!cRoot.CanHaveChild)
                cRoot = cRoot.Parent;
            ScenePresenter.SetCurrentRoot(cRoot);

            if (!cRoot.OpticalEntity.AddChild(newOSObject.OpticalEntity)) return false;

            cRoot.AddChild(newOSObject);
            return true;
        }

        internal void RemoveOpticalSceneObject(OpticalSceneObject oldOSObject)
        {
            if (oldOSObject.Root != SceneObjectTree) return;

            if (oldOSObject.Parent != null)
                ScenePresenter.SetCurrentRoot(oldOSObject.Parent);

            oldOSObject.Dispose();
        }
        private void AddRegistrators(List<OpticalSceneObject> regs)
        {
            foreach (var element in regs)
                if (!Registrators.Contains(element))
                    Registrators.Add(element);
        }
        private void RemoveRegistrators(List<OpticalSceneObject> regs)
        {
            foreach (var element in regs)
                Registrators.Remove(element);
        }
        private void AddSources(List<OpticalSceneObject> srcs)
        {
            foreach (var element in srcs)
                if (!Sources.Contains(element))
                    Sources.Add(element);
        }
        private void RemoveSources(List<OpticalSceneObject> srcs)
        {
            foreach (var element in srcs)
                Sources.Remove(element);
        }

        #endregion Methods

        #region Events

        protected void RaiseEvent(Object id, Object p)
        {
            var ocargs = new ObjectChangedEventArgs() { Source = this, ID = id, EventParams = p };

            if (ObjectChanged != null)
                ObjectChanged(this, ocargs);
        }

        private void OnOSObjectChanged(object source, RoutedChangeEventArgs args)
        {
            if (args.ID is OpticalSceneObject.ChildrenEventIDs)
            {
                var targetObj = args.EventParams as OpticalSceneObject;

                if (targetObj != null)
                {
                    var sources = targetObj.GetSources();
                    var registrators = targetObj.GetRegistrators();

                    switch ((OpticalSceneObject.ChildrenEventIDs) args.ID)
                    {
                        case OpticalSceneObject.ChildrenEventIDs.Add:
                            AddSources(sources);
                            AddRegistrators(registrators);

                            RaiseEvent(EventIDs.Add, targetObj);
                            break;
                        case OpticalSceneObject.ChildrenEventIDs.Remove:
                            RemoveSources(sources);
                            RemoveRegistrators(registrators);

                            RaiseEvent(EventIDs.Delete, targetObj);
                            break;
                    }
                }
            }
        }

        #endregion Events

        [OnDeserialized]
        internal void OnDeserializedMethod(StreamingContext context)
        {
            scenePresenter = new OpenGLOpticalScene(this);
        }
    }
}
