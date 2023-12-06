using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using Common;
using MathKernel.VectorAlgebra;
using ObjectModel.OpticalSceneEntities.OriginEntities;
using ObjectModel.SceneGraph;
using SharpGL;
using ViewModelCore.Common;

namespace ViewModelCore.OpticalSceneEntities.SceneView
{
    public class RenderSceneObject : IGLRenderable
    {
        private HashSet<RenderObject> renderObjects;

        public List<RenderSceneObject> Children { get; set; }
        internal RenderSceneObject Parent { get; set; }
        internal OpticalSceneObject OSObject { get; private set; }
        public RenderObject BoundingBox { get; set; }
        // ReSharper disable once MemberCanBeInternal
        public RenderOptionsS RenderOptions { get; set; }

        internal RenderSceneObject()
        {
            renderObjects = new HashSet<RenderObject>();
            Children = new List<RenderSceneObject>();
        }

        internal RenderSceneObject(OpticalSceneObject osObj)
            : this()
        {
            if (osObj == null) return;
            OSObject = osObj;
            WeakEventManager<OpticalSceneObject, ObjectChangedEventArgs>.AddHandler(OSObject, "ObjectChanged", OnOSObjectChanged);
        }

        internal void Hit(SortedList<double, RenderSceneObject> hittedObjects, OpenGL gl, int x, int y, MathKernel.VectorAlgebra.Vector ray)
        {
            //if (OSObject.OriginObject is ControlPoint || BoundingBox.Hits(gl, x, y, ray))
            {
                if (renderObjects.Count != 0)
                {
                    SortedList<double, RenderObject> hittedRenderObjects = new SortedList<double, RenderObject>(renderObjects.Count);

                    foreach (var obj in renderObjects)
                    {
                        try
                        {
                            obj.Hit(hittedRenderObjects, gl, x, y, ray);
                        }
                        catch { }
                    }

                    if (hittedRenderObjects.Count > 0)
                        try
                        {
                            hittedObjects.Add(hittedRenderObjects.First().Key, this);
                        }catch{}
                }
                foreach (var child in Children)
                    child.Hit(hittedObjects, gl, x, y, ray);
            }
        }

        internal void AddChildren(RenderSceneObject obj)
        {
            obj.Parent = this;
            Children.Add(obj);
        }

        internal void AddRenderObject(RenderObject obj)
        {
            renderObjects.Add(obj);
        }
        public void AddRenderObjects(List<RenderObject> objs)
        {
            foreach(var value in objs)
                AddRenderObject(value);
        }

        internal RenderSceneObject FindOSObject(OpticalSceneObject obj2Find)
        {
            if (OSObject == obj2Find) return this;

            foreach (var child in Children)
            {
                var curRes = child.FindOSObject(obj2Find);
                if (curRes != null) return curRes;
            }

            return null;
        }

        public void Render(OpenGL gl)
        {
            foreach (var value in renderObjects)
            {
                value.Render(gl);
            }
            foreach (var child in Children)
                child.Render(gl);
        }

        internal List<RenderSceneObject> UpdateLayout(LayoutRenderOptions layoutOpts)
        {
            foreach(var value in renderObjects)
                value.UpdateLayout(layoutOpts);
            return Children;
        }

        internal void UpdateLayout()
        {
            if (OSObject != null)
                OSObject.UpdateRenderLayout(this, RenderOptions);
            else
                foreach (var child in Children)
                    child.UpdateLayout();
        }

        internal void Dispose()
        {
            // remove listening
            if (OSObject != null)
                WeakEventManager<OpticalSceneObject, ObjectChangedEventArgs>.RemoveHandler(OSObject, "ObjectChanged",
                    OnOSObjectChanged);

            // Dispose render objects
            foreach (var obj in renderObjects)
                obj.Dispose();

            // Dispose children
            while (Children.Count > 0)
                Children[0].Dispose();

            // Remove self from parent
            try
            {
                Parent.Children.Remove(this);
            }
            catch { }
        }

        #region Events

        private void OnOSObjectChanged(object source, ObjectChangedEventArgs args)
        {
            if (Equals(args.ID, OpticalSceneObject.EventIDs.Layout))
                OSObject.UpdateRenderLayout(this, RenderOptions);
            else if (Equals(args.ID, OpticalSceneObject.EventIDs.Geometry) || Equals(args.ID, OpticalSceneObject.EventIDs.CS))
                Dispose();
            else if (OpticalSceneObject.ChildrenEventIDs.Remove.Equals(args.ID))
            {
                var osoChild = args.EventParams as OpticalSceneObject;
                if (osoChild != null)
                {
                    var rso = FindOSObject(osoChild);
                    if (rso != null) rso.Dispose();
                }
            }
        }

        #endregion Events
    }
}
