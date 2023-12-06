using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Runtime.Serialization;
using System.Text.RegularExpressions;
using System.Windows;
using Common;
using Common.Memento;
using Common.TypeExtensions;
using MathKernel;
using MathKernel.VectorAlgebra;
using Newtonsoft.Json;
using ObjectModel.Common;
using ObjectModel.OpticalSceneEntities;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ObjectModel.OpticalSceneEntities.OriginEntities.Surfaces;
using ObjectModel.Optimization.ControlPoints;
using ObjectModel.SceneGraph;
using ViewModelCore.Common;
using ViewModelCore.Controllers;
using ViewModelCore.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.OpticalSceneEntities.SceneView;
using ViewModelCore.Optimization.ControlPoints;

namespace ViewModelCore.OpticalSceneEntities
{
    [Serializable][JsonObject(MemberSerialization.Fields)]
    public class OpticalSceneObject : INotifyPropertyChanged, IInstanceOfList, INotifyDataErrorInfo
    {
        #region Static fields

        internal static int PointCounter = 1;
        internal static int SurfaceCounter = 1;
        internal static int SolidCounter = 1;
        internal static int GroupCounter = 1;
        internal static int SourceCounter = 1;
        internal static int RegistratorCounter = 1;

        #endregion Static fields

        #region Fields
        private string name = "Unknown";

        private readonly IDrawable originObject;
        private readonly OpticalEntity opticalEntity;
        private readonly List<OpticalSceneObject> childrenObjects = new List<OpticalSceneObject>();

        [field:NonSerialized, ReferenceMemorizableAttribute]
        private OpticalSceneObject parent;

        private bool isHidden = false;
        private bool isSelected = false;
        private bool isRoot = false;
        private bool hasControlPoints = false;
        private double tolerance = 0.01;

        [field: NonSerialized] private List<TriangleSet> surfaces;
        [field: NonSerialized] private List<Points> contour;

        [field: NonSerialized] private bool hasErrors;

        #endregion Fields

        [field: NonSerialized] public event EventHandler<ObjectChangedEventArgs> ObjectChanged;
        [field: NonSerialized] public event EventHandler<RoutedChangeEventArgs> RoutedChange;
        [field: NonSerialized] public event PropertyChangedEventHandler PropertyChanged;
        [field: NonSerialized] public event EventHandler<DataErrorsChangedEventArgs> ErrorsChanged;

        internal enum EventIDs { Geometry, CS, Layout, Disposed }
        internal enum ChildrenEventIDs { Add, Remove, Replace }

        #region Properties

        internal OpticalSceneObject Root
        {
            get
            {
                var res = this;
                while (res.Parent != null)
                    res = res.Parent;
                return res;
            }
        }

        internal OpticalSceneObject Parent
        {
            get { return parent; }
            private set { parent = value; }
        }

        public Object OriginObjectView
        {
            get
            {
                if (OpticalEntity != null)
                {
                    return OpticalEntityView.Wrap(OpticalEntity);
                }
                if (OriginObject is BezierControlPoint3D)
                {
                    return new BezierControlPoint3DView((BezierControlPoint3D)OriginObject, CSInGlobal);
                }

                throw new Exception(); // TODO
            }
        }

        public OpticalEntity OpticalEntity
        {
            get
            {
                return opticalEntity;
            }
        }

        public List<OpticalSceneObject> Children
        {
            get { return childrenObjects; }
        }

        internal CoordinateSystem CSInGlobal
        {
            get
            {
                var oe = OpticalEntity ?? Parent.OpticalEntity;
                if (oe != null) return oe.CSInGlobal;
                return CoordinateSystem.Global();
            }
        }

        public double Tolerance
        {
            get { return tolerance; }
            set
            {
                if (tolerance == value) return;

                SetTolerance(value);
                RaiseEvents(EventIDs.Geometry);
            }
        }

        public bool IsPoint
        {
            get { return OriginObject is ControlPoint; }
        }

        internal bool IsHidden
        {
            get { return isHidden; }
            set
            {
                if (value == isHidden) return;

                isHidden = value;
                if (value) isSelected = false;

                RaiseEvents(EventIDs.Layout);
                OnPropertyChanged();
            }
        }

        internal bool IsSelected
        {
            get { return isSelected; }
            set
            {
                if (value == isSelected) return;

                isSelected = value;
                RaiseEvents(EventIDs.Layout);
                OnPropertyChanged();
            }
        }

        internal bool IsRoot
        {
            get { return isRoot; }
            set
            {
                if (value == isRoot) return;

                isRoot = value;
                RaiseEvents(EventIDs.Layout);
            }
        }

        public bool IsRegistrator
        {
            get { return OriginObject is Registrator; }
        }

        public virtual bool IsSource
        {
            get
            {
                if (OpticalEntity == null) return false;

                return OpticalEntity.IsSource;
            }
        }

        public bool IsSurface
        {
            get
            {
                return OpticalEntity != null && OpticalEntity.IsSurface;
            }
        }

        public bool IsSolid
        {
            get { return OriginObject is SolidOpticalEntity; }
        }
        public bool IsGroup
        {
            get { return OpticalEntity != null && OpticalEntity.OEOrigin == null; }
        }

        internal bool CanHaveChild
        {
            get { return IsGroup && !IsSolid; }
        }

        internal bool IsRootAncestor
        {
            get
            {
                if (Children.Count == 0) return false;

                foreach (var child in Children)
                    if (child.IsRoot || child.IsRootAncestor) return true;

                return false;
            }
        }

        public string Name
        {
            get { return name; }
            set
            {
                if (Equals(name, value)) return;

                name = value;
                OnPropertyChanged();
            }
        }

        internal IDrawable OriginObject
        {
            get { return originObject; }
        }

        internal bool HasControlPoints
        {
            get { return hasControlPoints; }
            private set { hasControlPoints = value; }
        }

        private List<TriangleSet> Surfaces
        {
            get { return surfaces ?? (surfaces = OriginObject.GetSurfaces(Tolerance)); }
        }

        private List<Points> Contour
        {
            get { return contour ?? (contour = OriginObject.GetContourGeometry(Tolerance)); }
        }

        public IEnumerable GetErrors(string propertyName)
        {
            throw new NotImplementedException();
        }

        public bool HasErrors { get { return hasErrors; } protected set { hasErrors = value; } }

        #endregion Properties

        #region Methods
        protected OpticalSceneObject()
        {
        }

        [JsonConstructor]
        private OpticalSceneObject(IDrawable originObject)
        {
            if (originObject == null)
            {
                return;
            }
            this.originObject = originObject;
            opticalEntity = originObject as OpticalEntity;

            InitName();

            if (OpticalEntity != null)
                WeakEventManager<OpticalEntity, ObjectChangedEventArgs>.AddHandler(OpticalEntity, "ObjectChanged", OnOriginChanged);
        }

        public OpticalSceneObject(OpticalEntity opticalEntity)
            : this(opticalEntity as IDrawable)
        {
            if (opticalEntity != null)
            {
                childrenObjects = new List<OpticalSceneObject>(opticalEntity.Children.Count);
                foreach (var value in opticalEntity.Children)
                    AddChild(new OpticalSceneObject(value));

                if (opticalEntity is I3DPointEditable)
                {
                    foreach (var point in (opticalEntity as I3DPointEditable).GetControlPoints())
                        AddChild(new OpticalSceneObject(point));
                    HasControlPoints = true;
                }

                if (opticalEntity.OEOrigin is I3DPointEditable)
                {
                    var origin = (I3DPointEditable) opticalEntity.OEOrigin;
                    var points = origin.GetControlPoints();
                    foreach (var point in points)
                        AddChild(new OpticalSceneObject(point));
                    HasControlPoints = true;
                }
            }
        }

        #region Public

        internal void Move(MathKernel.VectorAlgebra.Point shiftVector, bool isAbs)
        {
            if (OpticalEntity != null)
                OpticalEntity.Move(shiftVector, isAbs);
            else if (OriginObject is BezierControlPoint3D)
            {
                var controlPoint = (BezierControlPoint3D)OriginObject;
                var oldPoint = controlPoint.Point;
                if (isAbs)
                    controlPoint.SetValues(shiftVector.X, shiftVector.Y, shiftVector.Z);
                else
                    controlPoint.SetValues(shiftVector.X + oldPoint.X, shiftVector.Y + oldPoint.Y, shiftVector.Z + oldPoint.Z);
            }
        }

        internal void Mirror(MathKernel.VectorAlgebra.Point normalStart, MathKernel.VectorAlgebra.Point normalEnd)
        {
            var norm = (normalEnd - normalStart).Normalize();

            if (OpticalEntity != null)
                OpticalEntity.Mirror(normalStart, norm);
        }

        internal void Rotate(MathKernel.VectorAlgebra.Point axisStart, MathKernel.VectorAlgebra.Point axisEnd, double angleRad)
        {
            var vec = (axisEnd - axisStart).Normalize();

            if (OpticalEntity != null)
                OpticalEntity.Rotate(axisStart, vec, angleRad);
        }

        internal OpticalSceneObject Copy()
        {
            var newName = Name;
            if (!Regex.IsMatch(newName, @"_copy$")) newName = newName + "_copy";

            if (OpticalEntity != null)
                return new OpticalSceneObject(OpticalEntity.Copy()) { Name = newName };

            throw new NotImplementedException();
        }

        public OpticalSceneObject Clone()
        {
            throw new NotImplementedException();
        }

        internal List<OpticalSceneObject> ToList()
        {
            List<OpticalSceneObject> result = new List<OpticalSceneObject>();
            FillList(result);
            return result;
        }

        internal List<OpticalSceneObject> GetSources()
        {
            List<OpticalSceneObject> result = new List<OpticalSceneObject>();
            FillSources(result);
            return result;
        }

        internal List<OpticalSceneObject> GetRegistrators()
        {
            List<OpticalSceneObject> result = new List<OpticalSceneObject>();
            FillRegistrators(result);
            return result;
        }

        internal static LayoutRenderOptions InitLayoutRenderOptions(RenderOptionsS renderOpts)
        {
            var res = new LayoutRenderOptions();

            // Point opts
            var pntOpts = new PointRenderOptions()
            {
                IsSmooth = renderOpts.RenderSmoothPoints,
                Size = (float)ResourceController.PointSize,
                Color = ColorExtensions.ToGLColor(ResourceController.PointColor)
            };
            if (renderOpts.RenderSelected)
            {
                pntOpts.Size = (float)ResourceController.SelectedPointSize;
                pntOpts.Color = ColorExtensions.ToGLColor(ResourceController.SelectedPointColor);
            }
            else if (renderOpts.RenderGhosted)
                pntOpts.Color = ColorExtensions.ToGLColor(ResourceController.DisabledColor);
            else if (renderOpts.RenderRoot)
                pntOpts.Color = ColorExtensions.ToGLColor(ResourceController.SceneRootColor);

            res.PointOptions = pntOpts;

            // Line opts
            var lineOpts = new LineRenderOptions()
            {
                IsSmooth = renderOpts.RenderSmoothPoints,
                Width = (float)ResourceController.LineWidth,
                Color = ColorExtensions.ToGLColor(ResourceController.LineColor)
            };
            if (renderOpts.RenderSelected)
            {
                lineOpts.Width = (float)ResourceController.SelectedLineWidth;
                lineOpts.Color = ColorExtensions.ToGLColor(ResourceController.SelectedLineColor);
            }
            else if (renderOpts.RenderGhosted)
                lineOpts.Color = ColorExtensions.ToGLColor(ResourceController.DisabledColor);
            else if (renderOpts.RenderRoot)
                lineOpts.Color = ColorExtensions.ToGLColor(ResourceController.SceneRootColor);

            res.LineOptions = lineOpts;

            // Line opts
            var contourOpts = new LineRenderOptions()
            {
                IsSmooth = renderOpts.RenderSmoothPoints,
                Width = (float)ResourceController.LineWidth,
                Color = ColorExtensions.ToGLColor(ResourceController.ContourLineColor)
            };

            res.ContourOptions = contourOpts;

            // Polygon opts
            pntOpts.Size = pntOpts.Size/2;
            lineOpts.Width = lineOpts.Width/3;
            var polyOpts = new PolygonRenderOptions()
            {
                FillColor = ColorExtensions.ToGLColor(ResourceController.SceneSurfaceColor),
                Mode = renderOpts.SurfaceRenderMode,
                PointOptions = pntOpts,
                LineOptions = lineOpts
            };
            if (renderOpts.RenderSelected)
                polyOpts.FillColor = ColorExtensions.ToGLColor(ResourceController.SceneSelectedSurfaceColor);
            else if (renderOpts.RenderGhosted)
                polyOpts.FillColor = ColorExtensions.ToGLColor(ResourceController.DisabledColor);
            else if (renderOpts.RenderRoot)
                polyOpts.FillColor = ColorExtensions.ToGLColor(ResourceController.SceneRootColor);

            res.PolygonOptions = polyOpts;

            return res;
        }

        #endregion Public

        #region Private

        protected void InitName()
        {
            if (IsPoint)
                Name = String.Format("Node {0}", PointCounter++);
            else if (IsSource)
                Name = String.Format("Source {0}", SourceCounter++);
            else if (IsRegistrator)
                Name = String.Format("Registrator {0}", RegistratorCounter++);
            else if (IsSolid)
                Name = String.Format("Solid {0}", SolidCounter++);
            else if (IsGroup)
                Name = String.Format("Group {0}", GroupCounter++);
            else if (IsSurface)
                Name = String.Format("Surface {0}", SurfaceCounter++);
        }

        private void FillList(List<OpticalSceneObject> objs)
        {
            objs.Add(this);
            for (int i = 0; i < Children.Count; i++)
                Children[i].FillList(objs);
        }

        private void FillSources(List<OpticalSceneObject> objs)
        {
            if (IsSource)
                objs.Add(this);
            else
                for (int i = 0; i < Children.Count; i++)
                    Children[i].FillSources(objs);
        }

        private void FillRegistrators(List<OpticalSceneObject> objs)
        {
            if (IsRegistrator)
                objs.Add(this);
            else
                for (int i = 0; i < Children.Count; i++)
                    Children[i].FillRegistrators(objs);
        }

        private RenderOptionsS FillRSObject(RenderSceneObject rsObject, RenderOptionsS renderOpts)
        {
            var curCoorSys = CSInGlobal;

            // Init layout options
            renderOpts = FitRenderOptions(renderOpts);
            var layoutOpts = InitLayoutRenderOptions(renderOpts);

            // Draw if is Control Point
            if (OriginObject is BezierControlPoint3D)
            {
                if (Parent.IsRoot)
                    rsObject.AddRenderObject(
                        new RenderPoint(curCoorSys.GetGlobalPoint((OriginObject as BezierControlPoint3D).CartesianPoint))
                        {
                            Options = layoutOpts.PointOptions,
                            SourceObject = OriginObject,
                            SceneObject = this
                        }
                        );
                else rsObject.Dispose();
                return renderOpts;
            }

            // Draw if is not Control Point
            // prepare surfaces
            if (Surfaces != null)
                foreach (var value in Surfaces)
                {
                    var globalVal = value.DeepClone() as TriangleSet;
                    globalVal.ToGlobal(curCoorSys);
                    rsObject.AddRenderObject(new RenderTriangles(globalVal)
                    {
                        Options = layoutOpts.PolygonOptions,
                        SourceObject = value,
                        SceneObject = this
                    }
                        );
                }

            // prepare contours
            if (renderOpts.RenderContours && Contour != null)
            {
                foreach (var value in Contour)
                {
                    var globalVal = curCoorSys.GetGlobalPoint(value);
                    rsObject.AddRenderObject(new RenderLine(globalVal)
                    {
                        Options = layoutOpts.ContourOptions,
                        IsContour = true,
                        SourceObject = value,
                        SceneObject = this
                    }
                        );
                }
            }

            return renderOpts;
        }

        private void SetTolerance(double value)
        {
            tolerance = value;
            RefreshGeometry();
            foreach (var child in childrenObjects)
                child.SetTolerance(value);
        }

        private void RefreshGeometry()
        {
            surfaces = null;
            contour = null;
        }

        #endregion Private

        #region Internal

        internal bool AddChild(OpticalSceneObject child)
        {
            if (childrenObjects.Contains(child)) return false;

            child.Parent = this;
            childrenObjects.Add(child);
            SubscribeToChild(child);

            if (OpticalEntity != null && child.OpticalEntity != null && !OpticalEntity.Children.Contains(child.OriginObject))
                OpticalEntity.AddChild(child.OpticalEntity);

            RaiseEvents(ChildrenEventIDs.Add, child);
            return true;
        }
        internal bool RemoveChild(OpticalSceneObject child)
        {
            if (!childrenObjects.Contains(child) || !childrenObjects.Remove(child)) return false;

            UnsubscribeChild(child);
            child.Parent = null;

            if (OpticalEntity != null && OpticalEntity.Children.Contains(child.OriginObject))
                OpticalEntity.RemoveChild(child.OpticalEntity);

            RaiseEvents(ChildrenEventIDs.Remove, child);
            return true;
        }

        internal void PrepareRender(RenderSceneObject parentRSObj, RenderOptionsS renderOpts)
        {
            if (IsHidden) return;

            RenderSceneObject currentRSObj = new RenderSceneObject(this) {RenderOptions = renderOpts};
            parentRSObj.AddChildren(currentRSObj);

            renderOpts = FillRSObject(currentRSObj, renderOpts);

            foreach (var child in childrenObjects)
                child.PrepareRender(currentRSObj, renderOpts);
        }
        internal void UpdateRenderFull(RenderSceneObject correspObj, RenderOptionsS renderOpts)
        {
            if (IsHidden) correspObj.Dispose();

            renderOpts = FillRSObject(correspObj, renderOpts);

            // Render or update children
            int rsCNum = correspObj.Children.Count;
            Dictionary<RenderSceneObject, bool> childrenBools = new Dictionary<RenderSceneObject, bool>(rsCNum);
            for (var i = 0; i < rsCNum; i++)
                childrenBools.Add(correspObj.Children[i], false);

            foreach (var child in Children)
            {
                var foundRSChild = false;

                for (var i = 0; i < rsCNum; i++)
                {
                    var rsChild = correspObj.Children[i];
                    if (rsChild.OSObject == child)
                    {
                        child.UpdateRenderFull(rsChild, renderOpts);
                        childrenBools[rsChild] = true;
                        foundRSChild = true;
                        break;
                    }
                }

                if (!foundRSChild)
                    child.PrepareRender(correspObj, renderOpts);
            }

            // Remove redundant render scene objects
            foreach (var value in childrenBools.Where(value => value.Value == false))
                value.Key.Dispose();

        }
        internal void UpdateRenderLayout(RenderSceneObject correspObj, RenderOptionsS renderOpts)
        {
            if (IsHidden)
            {
                correspObj.Dispose();
                return;
            }

            // Init layout options
            renderOpts = FitRenderOptions(renderOpts);
            var layoutOpts = InitLayoutRenderOptions(renderOpts);
            var renderedChildren = correspObj.UpdateLayout(layoutOpts);
            foreach (var child in renderedChildren)
            {
                int pos = childrenObjects.IndexOf(child.OSObject);
                if (pos == -1) child.Dispose();
                else childrenObjects[pos].UpdateRenderLayout(child, renderOpts);
            }
        }
        internal RenderOptionsS FitRenderOptions(RenderOptionsS obj2Fit)
        {
            if (IsSelected)
            {
                obj2Fit.RenderSelected = true;
            }
            if (IsRoot)
            {
                obj2Fit.RenderGhosted = false;
                obj2Fit.RenderRoot = true;
            }
            else
                obj2Fit.RenderRoot = false;

            return obj2Fit;
        }

        internal void Dispose()
        {
            Parent.RemoveChild(this);

            RaiseEvents(EventIDs.Disposed);
        }

        internal static void ResetCounters()
        {
            PointCounter = 1;
            SurfaceCounter = 1;
            SolidCounter = 1;
            GroupCounter = 1;
            SourceCounter = 1;
            RegistratorCounter = 1;
        }

        #endregion Internal

        #endregion Methods

        #region Events

        protected void RaiseEvents(Object id)
        {
            RaiseEvents(id, null);
        }

        protected void RaiseEvents(Object id, Object Params)
        {
            var ocargs = new ObjectChangedEventArgs() { Source = this, ID = id, EventParams = Params };
            var rcargs = new RoutedChangeEventArgs() { Source = this, OriginalSource = this, ID = id, Handled = false, EventParams = Params };

            if (ObjectChanged != null)
                ObjectChanged(this, ocargs);
            RaiseRoutedChangeEvent(rcargs);
        }

        private void OnOriginChanged(object source, ObjectChangedEventArgs args)
        {
            if (Equals(args.ID, OpticalEntity.EventIDs.Geometry))
            {
                RefreshGeometry();
                RaiseEvents(EventIDs.Geometry);
            }
            else if (Equals(args.ID, OpticalEntity.EventIDs.CS))
                RaiseEvents(EventIDs.CS);
        }

        protected virtual void OnPropertyChanged([CallerMemberName] string propertyName = null)
        {
            var handler = PropertyChanged;
            if (handler != null) handler(this, new PropertyChangedEventArgs(propertyName));
        }

        // ROUTED EVENTS

        protected void SubscribeToChildren(List<OpticalSceneObject> Children)
        {
            foreach (var value in Children)
                SubscribeToChild(value);
        }

        protected void SubscribeToChild(OpticalSceneObject Child)
        {
            WeakEventManager<OpticalSceneObject, RoutedChangeEventArgs>.AddHandler(Child, "RoutedChange",
                OnObjectChangedInternal);
        }

        protected void UnsubscribeChild(OpticalSceneObject Child)
        {
            WeakEventManager<OpticalSceneObject, RoutedChangeEventArgs>.RemoveHandler(Child, "RoutedChange",
                OnObjectChangedInternal);
        }

        protected void RaiseRoutedChangeEvent(RoutedChangeEventArgs args)
        {
            if (RoutedChange != null)
                RoutedChange(this, args);
        }

        private void OnObjectChangedInternal(object source, RoutedChangeEventArgs args)
        {
            OnObjectChanged(source, args);

            if (args.Handled) return;
            RaiseRoutedChangeEvent(args);
        }

        protected virtual void OnObjectChanged(object source, RoutedChangeEventArgs args)
        {
        }

        #endregion Events

        [OnDeserialized]
        internal void OnDeserializedMethod(StreamingContext context)
        {
            foreach (var entity in childrenObjects)
            {
                entity.Parent = this;
            }
        }
    }
}
