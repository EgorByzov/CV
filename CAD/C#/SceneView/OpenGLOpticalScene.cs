using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using WPoint = System.Windows.Point;
using System.Windows.Input;
using Common;
using MathKernel;
using Newtonsoft.Json;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ObjectModel.SceneGraph;
using SharpGL;
using ViewModelCore.Common;
using ObjectModel.OpticalSceneEntities;
using ViewModelCore.Controllers;
using Common.TypeExtensions;
namespace ViewModelCore.OpticalSceneEntities.SceneView
{
    public class OpenGLOpticalScene : OpenGLRenderScene
    {
        private Dictionary<OpticalEntity, CoordinateSystem> _globalCSTableCash;
        
        public event EventHandler<ObjectChangedEventArgs> ObjectChanged;
        public enum EventIDs { Selection, Root };

        public Object SelectedObject
        {
            get
            {
                if (SelectedObjects.Count == 0)
                    return null;
                if (SelectedObjects.Count == 1)
                    return SelectedObjects.First();
                return SelectedObjects;
            }
        }

        internal HashSet<OpticalSceneObject> SelectedObjects { get; private set; }
        internal HashSet<OpticalSceneObject> HiddenObjects { get; private set; }
        internal OpticalSceneObject CurrentRoot { get; private set; }
        public OpticalScene OS { get; private set; }

        public OpenGLOpticalScene() :this(new OpticalScene()){}

        [JsonConstructor]
        public OpenGLOpticalScene(OpticalScene os)
        {
            OS = os;
            SelectedObjects = new HashSet<OpticalSceneObject>();
            HiddenObjects = new HashSet<OpticalSceneObject>();
            SetCurrentRoot(OS.SceneObjectTree);

            SubscribeTree(OS.SceneObjectTree);
            OS.SceneObjectTree.RoutedChange += OnOSObjectChanged;

            InitColorPalette();
            ApplicationPreferences.StaticPropertyChanged += (o, e) =>
            {
                if (e.PropertyName.Equals("ColorPalette"))
                    InitColorPalette();
            };

            Render();
        }

        internal void SetCurrentRoot(OpticalSceneObject obj)
        {
            if (CurrentRoot == obj || obj == null) return;

            bool IsRenderNeeded = obj.HasControlPoints;

            if (CurrentRoot != null)
            {
                if (obj.Children == null || obj.Children.Count == 0) return;

                CurrentRoot.IsRoot = false;
                IsRenderNeeded = IsRenderNeeded | CurrentRoot.HasControlPoints;
            }

            UnselectAll();
            CurrentRoot = obj;
            obj.IsRoot = true;

            RaiseEvent(EventIDs.Root, obj);
            if (IsRenderNeeded) RenderObject(CurrentRoot);
        }
        
        internal bool IsMoveAvailable()
        {
            if (SelectedObjects.Count < 1) return false;
            return true;
        }
        internal void MoveSelectedObjects(ImmutablePointView shiftVector, bool isAbs)
        {
            if (!IsMoveAvailable()) return;

            foreach (var el in SelectedObjects)
                el.Move(shiftVector.ModelObject, isAbs);
        }
        internal bool IsMirrorAvailable()
        {
            if (SelectedObjects.Count < 1) return false;
            return true;
        }
        internal void MirrorSelectedObjects(ImmutablePointView normalStart, ImmutablePointView normalEnd)
        {
            if (!IsMirrorAvailable()) return;

            foreach (var el in SelectedObjects)
                el.Mirror(normalStart.ModelObject, normalEnd.ModelObject);
        }
        internal bool IsRotateAvailable()
        {
            if (SelectedObjects.Count < 1) return false;
            return true;
        }
        internal void RotateSelectedObjects(ImmutablePointView axisStart, ImmutablePointView axisEnd, double angleRad)
        {
            if (!IsRotateAvailable()) return;

            foreach (var el in SelectedObjects)
                el.Rotate(axisStart.ModelObject, axisEnd.ModelObject, angleRad);
        }
        internal bool IsDuplicateAvailable()
        {
            if (SelectedObjects.Count < 1) return false;
            foreach (var el in SelectedObjects)
                if (el.OpticalEntity == null) return false;
            return true;
        }
        internal void DuplicateSelectedObjects()
        {
            if (!IsDuplicateAvailable()) return;

            foreach (var el in SelectedObjects)
                OS.AddOpticalSceneObject(el.Copy());
        }

        internal bool IsDeleteAvailable()
        {
            if (SelectedObject == null) return false;

            return !CurrentRoot.IsSolid && !CurrentRoot.IsSurface;
        }
        internal void DeleteSelectedObject()
        {
            if (!IsDeleteAvailable()) return;

            foreach (var el in SelectedObjects)
                el.Dispose();

            UnselectAll();
        }
        internal bool IsGroupAvailable()
        {
            return IsDeleteAvailable() && SelectedObjects.Count > 1;
        }
        internal void GroupSelectedObjects()
        {
            if (!IsGroupAvailable()) return;

            var selObjs = new List<OpticalSceneObject>(SelectedObjects);
            DeleteSelectedObject();

            var groupEntity = new OpticalEntity();
            var groupObj = new OpticalSceneObject(groupEntity);
            for (int i = 0; i < selObjs.Count; i++)
                groupObj.AddChild(selObjs[i]);

            OS.AddOpticalSceneObject(groupObj);
            ExclusiveSelect(groupObj);
        }
        internal bool IsUnGroupAvailable()
        {
            if (SelectedObjects.Count == 0) return false;

            foreach (var el in SelectedObjects)
                if (!el.IsGroup || el.IsSurface || el.IsSolid || el.Children.Count == 0) return false;

            return true;
        }
        internal void UnGroupSelectedObjects()
        {
            if (!IsUnGroupAvailable()) return;

            var selObjs = new List<OpticalSceneObject>(SelectedObjects);
            DeleteSelectedObject();

            foreach (var elOSO in selObjs)
                UnGroupObject(elOSO);
        }
        internal bool IsExplodeAvailable()
        {
            if (SelectedObjects.Count != 1) return false;
            return SelectedObjects.First().IsSolid;
        }
        internal void ExplodeSelectedObject()
        {
            if (!IsExplodeAvailable()) return;

            var solid = SelectedObjects.First();
            DeleteSelectedObject();

            UnGroupObject(solid);
        }
        private void UnGroupObject(OpticalSceneObject osObj)
        {
            var orig = osObj.OpticalEntity;
            if (orig != null)
            {
                var baseCS = orig.CS;

                foreach (var el in osObj.Children)
                {
                    var elOrig = el.OpticalEntity;
                    if (elOrig != null)
                        elOrig.CS = baseCS.GetGlobalCS(elOrig.CS);
                }
            }

            foreach (var el in osObj.Children)
                OS.AddOpticalSceneObject(el);
        }
        
        private void Render()
        {
            RenderSceneObject rs = new RenderSceneObject();
            var rsOpts = RenderOptionsS.Default;
            rsOpts.SurfaceRenderMode = OpenGL.GL_FILL;

            OS.SceneObjectTree.PrepareRender(rs, rsOpts);

            Render(rs);
        }

        private void RenderObject(OpticalSceneObject objToRender)
        {
            if (renderTree == null) return;

            RenderSceneObject rObj = renderTree.FindOSObject(objToRender);
            if (rObj == null)
            {
                rObj = renderTree.FindOSObject(objToRender.Parent);
                if (rObj == null) throw new Exception("This object is not yet in the tree");

                objToRender.PrepareRender(rObj, objToRender.Parent.FitRenderOptions(rObj.RenderOptions));
            }
            else
            {
                if (rObj.Parent.OSObject != objToRender.Parent) throw new Exception("This object is allready rendered");

                objToRender.UpdateRenderFull(rObj, rObj.RenderOptions);
            }
        }

        public void ZoomSelected()
        {
            ZoomAsync(SelectedObjects.ToList());
        }
        public void ZoomAll()
        {
            ZoomAsync(CurrentRoot.Children.Where(o => !o.IsRegistrator).ToList());
        }

        private async void ZoomAsync(List<OpticalSceneObject> objects)
        {
            if (objects.Count == 0) return;

            var border = GeometryBorders.Default;

            for (int i = 0; i < objects.Count; i++)
            {
                var curBorder = objects[i].OriginObject.Borders;
                var curOE = objects[i].OpticalEntity;
                if (curOE != null)
                    curBorder = curBorder.GetGlobal(curOE.CSInGlobal);
                border.UpdateBorders(curBorder);
            }

            await Task.Run(() => ZoomInMotion(border));
        }

        public async void SetViewAsync(StandardViews type)
        {
            await Task.Run(() => SetViewInMotion(type));
        }

        public void LeftBtnClickHandler(WPoint pos)
        {
            OpticalSceneObject selectableParent;
            var hittedObjects = HitScene(pos);

            ////// Click on empty space
            if (hittedObjects.Count == 0)
            {

                if ((Keyboard.Modifiers | ModifierKeys.None) == ModifierKeys.None)
                    UnselectAll();
                return;
            }

            ////// Click on object
            var firstHittedObject = hittedObjects.First().Value.OSObject;
            if (firstHittedObject == CurrentRoot) return;

            // get sellectable parent
            selectableParent = GetSelectableParent(firstHittedObject);
            if (selectableParent == null) return;

            // Single click
            if (Keyboard.IsKeyDown(Key.LeftShift) || Keyboard.IsKeyDown(Key.RightShift))
            {
                ExtraSelect(selectableParent);
                return;
            }

            if (Keyboard.IsKeyDown(Key.LeftCtrl) || Keyboard.IsKeyDown(Key.RightCtrl))
            {
                ExpulsiveSelect(selectableParent);
                return;
            }

            ExclusiveSelect(selectableParent);
        }

        public void LeftBtnDblClickHandler(WPoint pos)
        {
            var hittedObjects = HitScene(pos);

            if (hittedObjects.Count == 0)
            {
                if (CurrentRoot.Parent == null) return;
                SetCurrentRoot(CurrentRoot.Parent);
                return;
            }

            ////// Click on object
            var firstHittedObject = hittedObjects.First().Value.OSObject;
            if (firstHittedObject == CurrentRoot) return;

            // get sellectable parent
            var selectableParent = GetSelectableParent(firstHittedObject);
            if (selectableParent == null) return;

            SetCurrentRoot(selectableParent);
        }

        public override SortedList<double, RenderSceneObject> HitScene(WPoint pos)
        {
            var currentRoot = renderTree.FindOSObject(CurrentRoot);
            return HitScene(pos, currentRoot);
        }

        public void KeyDownHandler(Key key)
        {
            switch (key)
            {
                case Key.W:
                    Zoom(1);
                    break;
                case Key.S:
                    Zoom(-1);
                    break;
                case Key.Up:
                    KeyRotate(-1, 0);
                    break;
                case Key.Down:
                    KeyRotate(+1, 0);
                    break;
                case Key.Left:
                    KeyRotate(0, -1);
                    break;
                case Key.Right:
                    KeyRotate(0, +1);
                    break;
                case Key.NumPad8:
                    ChangeFieldOfView(1);
                    break;
                case Key.NumPad2:
                    ChangeFieldOfView(-1);
                    break;
            }
        }

        public void UnhideAll()
        {
            foreach (var obj in HiddenObjects)
            {
                obj.IsHidden = false;
            }
            HiddenObjects.Clear();
        }
        public void HideSelected()
        {
            foreach (var obj in SelectedObjects)
            {
                obj.IsHidden = true;
                HiddenObjects.Add(obj);
            }
            UnselectAll();
        }
        // The only one selected object
        internal bool ExclusiveSelect(OpticalSceneObject selectedObject)
        {
            UnselectAll(true);
            SelectObject(selectedObject);
            RaiseEvent(EventIDs.Selection, selectedObject);

            return true;
        }
        // Plus selected object
        internal bool ExtraSelect(OpticalSceneObject selectedObject)
        {
            SelectObject(selectedObject);
            RaiseEvent(EventIDs.Selection, selectedObject);

            return true;
        }
        // Minus selected object
        internal bool ExpulsiveSelect(OpticalSceneObject selectedObject)
        {
            selectedObject.IsSelected = false;
            SelectedObjects.Remove(selectedObject);
            RaiseEvent(EventIDs.Selection, selectedObject);

            return true;
        }
        //
        private void SelectObject(OpticalSceneObject selectedObject)
        {
            if (!selectedObject.Parent.IsRoot)
                SetCurrentRoot(selectedObject.Parent);

            selectedObject.IsSelected = true;
            SelectedObjects.Add(selectedObject);
        }
        // Clear selection
        private void UnselectAll(bool isSilentMode = false)
        {
            foreach (var obj in SelectedObjects)
            {
                obj.IsSelected = false;
            }
            SelectedObjects.Clear();

            if (!isSilentMode)
                RaiseEvent(EventIDs.Selection, null);
        }

        private OpticalSceneObject GetSelectableParent(OpticalSceneObject selectedObject)
        {
            var iter = selectedObject;
            while (iter != null && iter.Parent != CurrentRoot)
            {
                iter = iter.Parent;
            }

            return iter;
        }

        private void SubscribeTree(OpticalSceneObject node)
        {
            SubscribeChild(node);
            foreach (var opticalSceneObject in node.Children)
                SubscribeTree(opticalSceneObject);
        }

        private void SubscribeChild(OpticalSceneObject child)
        {
            WeakEventManager<OpticalSceneObject, ObjectChangedEventArgs>.AddHandler(child, "ObjectChanged", (o, e) =>
            {
                if (Equals(e.ID, OpticalSceneObject.EventIDs.Geometry))
                    RenderObject((OpticalSceneObject)o);
            });
        }

        protected void RaiseEvent(Object id, Object p)
        {
            var ocargs = new ObjectChangedEventArgs() { Source = this, ID = id, EventParams = p };

            if (ObjectChanged != null)
                ObjectChanged(this, ocargs);
        }

        private void OnOSObjectChanged(object source, RoutedChangeEventArgs args)
        {
            if (args.ID is OpticalSceneObject.EventIDs)
            {
                if (Equals(args.ID, OpticalSceneObject.EventIDs.Geometry))
                    RenderObject(args.Source as OpticalSceneObject);
                if (Equals(args.ID, OpticalSceneObject.EventIDs.CS))
                {
                    RenderObject(args.Source as OpticalSceneObject);
                }
            }
            else if (args.ID is OpticalSceneObject.ChildrenEventIDs)
            {
                var targetObj = args.EventParams as OpticalSceneObject;

                if (targetObj != null)
                {
                    switch ((OpticalSceneObject.ChildrenEventIDs) args.ID)
                    {
                        case OpticalSceneObject.ChildrenEventIDs.Add:                            
                            SubscribeChild(targetObj);
                            RenderObject(targetObj);
                            break;
                    }
                }
            }
        }

        private void InitColorPalette()
        {
            BackgroundColor = ColorExtensions.ToGLColor(ResourceController.SceneBackgroundColor);
            AxisXColor = ColorExtensions.ToGLColor(ResourceController.SceneAxisXColor);
            AxisYColor = ColorExtensions.ToGLColor(ResourceController.SceneAxisYColor);
            AxisZColor = ColorExtensions.ToGLColor(ResourceController.SceneAxisZColor);
            GridColor = ColorExtensions.ToGLColor(ResourceController.SceneNetColor);
            TargetColor = ColorExtensions.ToGLColor(ResourceController.SceneSelectedSurfaceColor); 

            if (renderTree != null)
                renderTree.UpdateLayout();
        }

    }
}
