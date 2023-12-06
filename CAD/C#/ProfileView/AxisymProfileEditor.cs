using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.ComponentModel;
using System.Linq;
using System.Runtime.CompilerServices;
using Common;
using MathKernel;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ObjectModel.OpticalSceneEntities.OriginEntities.Lines;
using ObjectModel.OpticalSceneEntities.OriginEntities.Surfaces;
using OxyPlot;
using OxyPlot.Wpf;
using ViewModelCore.Controllers;
using ViewModelCore.Optimization.ControlPoints;

namespace ViewModelCore.OpticalSceneEntities.ProfileView
{
    public enum ProfileEditorMode
    {
        Select,
        AddNode,
        RemoveNode,
        OptimizeNode,
        LockNode,
        Spline2Segline,
        Segline2Spline,
        CutLine,
        View
    };

    public class AxisymProfileEditor : INotifyPropertyChanged
    {
        # region Fields

        // static
        private static double HitTestTolerance = 5;
        private static double BindingThreshold = 7;

        // flags
        private bool _isSpeciallyEditable;
        private bool _isInteractiveMove;
        private bool _isEndPointBindingOn;
        private bool _isNearPointBindingOn;
        private bool _isControlPointBindingOn;
        private ProfileEditorMode _editorMode;

        // common line series
        private LineSeries _selectedPointMarker;

        // main objects
        private Plot _plotView;
        private OpticalSceneObject _editedObject;

        // after object processed fields
        private List<Surface> _passiveSurfaces;
        private List<AxisymLine2DSurface> _activeSurfaces;
        private Dictionary<Line2D, AxisymLine2DSurface> _lineToSurfaceMap;
        private Dictionary<Line2D, Line2DOxyController> _lineToControllerMap;
        private Dictionary<Surface, Line2DOxyController> _passiveSurfaceToControllerMap;
        private Dictionary<Surface, CoordinateSystem> _surfaceToCSMap;
        private List<AxisymSurfaceBinder> _activeSurfacesLinks;

        // editor mechanics
        private Line2DOxyController _selectedLine;
        private ControlPoint2DOxyController _selectedPoint;

        private bool _isStartingPointMove;
        private bool _isOnPointMove;
        private ScreenPoint _targetPoint;
        private List<LineSeries> _pointMoveAnimations;

        # endregion Fields

        # region Properties

        public OpticalSceneObject SceneObject { get; private set; }

        public bool IsSolid
        {
            get { return (_editedObject != null) && _editedObject.IsSolid; }
        }

        public bool IsSpeciallyEditable { get { return _isSpeciallyEditable; } private set { _isSpeciallyEditable = value; } }

        public bool IsInteractiveMove { get { return _isInteractiveMove; } set { _isInteractiveMove = value; } }

        public bool HasEditableLines { get { return _activeSurfaces.Count > 0; } }

        // Point Binding On Point Move
        public bool IsPointBindingOn { get { return IsEndPointBindingOn || IsNearPointBindingOn || IsControlPointBindingOn; } }
        public bool IsEndPointBindingOn { get { return _isEndPointBindingOn; } set { _isEndPointBindingOn = value; } }
        public bool IsNearPointBindingOn { get { return _isNearPointBindingOn; } set { _isNearPointBindingOn = value; } }
        public bool IsControlPointBindingOn { get { return _isControlPointBindingOn; } set { _isControlPointBindingOn = value; } }

        public ProfileEditorMode EditorMode {
            get { return _editorMode; }
            set
            {
                if (EditorMode == value) return;

                StopActions();
                _editorMode = value;
                foreach (var controller in _lineToControllerMap)
                    controller.Value.DisplayMode = value;
            }
        }

        public ControlPoint2DView SelectedPointView { get; private set; }

        internal Line2DOxyController SelectedLine
        {
            get { return _selectedLine; }
            set
            {
                if (SelectedLine == value) return;

                if (SelectedLine != null) SelectedLine.IsSelected = false;
                if (SelectedControlPoint != null && SelectedControlPoint.ParentLineController != value) SelectedControlPoint = null;
                _selectedLine = value;
                if (SelectedLine != null) SelectedLine.IsSelected = true;
                OnPropertyChanged();
            }
        }

        internal ControlPoint2DOxyController SelectedControlPoint
        {
            get { return _selectedPoint; }
            set
            {
                if (SelectedControlPoint == value) return;

                if (SelectedControlPoint != null) SelectedControlPoint.IsSelected = false;
                _selectedPoint = value;
                if (SelectedControlPoint != null) SelectedControlPoint.IsSelected = true;

                if (SelectedControlPoint == null) SelectedPointView = null;
                else SelectedPointView = ControlPoint2DView.Wrap(SelectedControlPoint.ControlPoint, SelectedControlPoint.CastCS);

                OnPropertyChanged();
            }
        }

        # endregion Properties

        private AxisymProfileEditor()
        {
            _isSpeciallyEditable = false;
            _isInteractiveMove = true;
            _passiveSurfaces = new List<Surface>();
            _activeSurfaces = new List<AxisymLine2DSurface>();
            _lineToSurfaceMap = new Dictionary<Line2D, AxisymLine2DSurface>();
            _lineToControllerMap = new Dictionary<Line2D, Line2DOxyController>();
            _passiveSurfaceToControllerMap = new Dictionary<Surface, Line2DOxyController>();
            _surfaceToCSMap = new Dictionary<Surface, CoordinateSystem>();
            _activeSurfacesLinks = new List<AxisymSurfaceBinder>();
            _pointMoveAnimations = new List<LineSeries>()
            {
                new LineSeries(){LineStyle = LineStyle.Dash, StrokeThickness = ResourceController.LineWidth, Color = ResourceController.DisabledColor},
                new LineSeries(){LineStyle = LineStyle.None, CanTrackerInterpolatePoints = false, MarkerType = MarkerType.Circle, MarkerStroke = ResourceController.ContourLineColor, MarkerSize = ResourceController.SelectedPointSize}
            };
            EditorMode = ProfileEditorMode.View;
        }

        public AxisymProfileEditor(OpticalSceneObject obj, Plot plot) : this()
        {
            //plot.Series.Clear();
            _plotView = plot;
            _editedObject = obj;
            if (obj != null)
                InitSceneObject();

            _plotView.Series.Add(_pointMoveAnimations[0]);
            _plotView.Series.Add(_pointMoveAnimations[1]);
        }

        public void NextCP()
        {
        }

        public void PreviousCP()
        {
        }

        public void HitEventHandler(HitTestArguments args)
        {
            if (EditorMode == ProfileEditorMode.View) return;
            _targetPoint = args.Point;

            double closestPointDist = Double.PositiveInfinity;
            double closestLineDist = Double.PositiveInfinity;
            ControlPoint2DOxyController closestPoint = null;
            Line2DOxyController closestLine = null;

            foreach (var element in _lineToControllerMap)
            {
                var curHitResults = element.Value.HitLineTest(args);
                if (curHitResults.IsLineHitted && curHitResults.Distance2Line() < closestLineDist)
                {
                    closestLine = curHitResults.HittedLine;
                    closestLineDist = curHitResults.Distance2Line();
                }
                if (curHitResults.IsControlPointHitted && curHitResults.Distance2ControlPoint() < closestPointDist)
                {
                    closestPoint = curHitResults.HittedControlPoint;
                    closestPointDist = curHitResults.Distance2ControlPoint();
                }
            }

            // Switch Editor Mode
            if (EditorMode == ProfileEditorMode.Select) // SELECT
            {
                if (closestPoint != null)
                {
                    SelectedControlPoint = closestPoint;
                    _isStartingPointMove = true;
                }
                else if (closestLine != null)
                    SelectedLine = closestLine;

                return;
            }
            if (EditorMode == ProfileEditorMode.OptimizeNode && closestPoint != null && !closestPoint.ControlPoint.IsLocked) // OPTIMIZE
            {
                closestPoint.ControlPoint.IsOptimizationOn = !closestPoint.ControlPoint.IsOptimizationOn;
                return;
            }
            if (EditorMode == ProfileEditorMode.LockNode && closestPoint != null && closestPoint.ControlPoint.IsUnlockAllowed) // LOCK
            {
                closestPoint.ControlPoint.IsLocked = !closestPoint.ControlPoint.IsLocked;
                return;
            }
            if (EditorMode == ProfileEditorMode.AddNode || EditorMode == ProfileEditorMode.RemoveNode)
            {
                if (closestLine != null && closestLineDist > closestPointDist)
                {
                    SelectedLine = closestLine;
                    return;
                }

                if (EditorMode == ProfileEditorMode.RemoveNode &&
                    closestPoint != null && closestPoint.ParentLineController == SelectedLine)
                {
                    SelectedLine.RemoveControlPoint(closestPoint);
                    return;
                }

                if (EditorMode == ProfileEditorMode.AddNode && SelectedLine != null && closestPoint == null)
                {
                    SelectedLine.AddControlPoint(Screen2PlotCoords(_targetPoint));
                    return;
                }
            }

        }
        public void OnPlotMouseUp()
        {
            if (_isOnPointMove)
                MovePoint(_targetPoint);
            StopActions();
        }
        public void OnPlotMouseMove(ScreenPoint mousePosition)
        {
            _targetPoint = mousePosition;

            if (_isStartingPointMove)
            {
                var cpPos = PlotCoords2ScenePoint(SelectedControlPoint.OxyPoint);
                if ((_targetPoint - cpPos).Length > HitTestTolerance)
                {
                    _isStartingPointMove = false;
                    _isOnPointMove = true;
                }
            }
            if (_isOnPointMove)
            {
                if (IsPointBindingOn)
                {
                    var dataPoint = Screen2PlotCoords(_targetPoint);
                    var closestPoints = new Dictionary<DataPoint, double>();

                    foreach (var pair in _lineToControllerMap)
                    {
                        var curLineController = pair.Value;
                        FillWithClosestPoints(dataPoint, closestPoints, curLineController);
                    }

                    foreach (var pair in _passiveSurfaceToControllerMap)
                    {
                        var curLineController = pair.Value;
                        FillWithClosestPoints(dataPoint, closestPoints, curLineController);
                    }
                    if (closestPoints.Count > 0)
                    {
                        var closestPoint = closestPoints.Aggregate((a, b) => b.Value < a.Value ? b : a);
                        var closestScreenPoint = PlotCoords2ScenePoint(closestPoint.Key);
                        if ((closestScreenPoint - _targetPoint).Length <= BindingThreshold)
                            _targetPoint = closestScreenPoint;
                    }
                }

                if (IsInteractiveMove)
                {
                    MovePoint(_targetPoint);
                    return;
                }

                var start = SelectedControlPoint.OxyPoint;
                var end = Screen2PlotCoords(_targetPoint);
                _pointMoveAnimations[0].ItemsSource = new List<DataPoint>()
                {
                    start, end
                };
                _pointMoveAnimations[1].ItemsSource = new List<DataPoint>()
                {
                    end
                };
            }
        }

        public void StopActions()
        {
            if (_isOnPointMove)
            {
                foreach (var line in _pointMoveAnimations)
                    line.ItemsSource = null;
            }

            _isStartingPointMove = false;
            _isOnPointMove = false;
        }
        public void UnselectAll()
        {
            SelectedLine = null;
            SelectedControlPoint = null;
        }

        private void MovePoint(ScreenPoint targetPoint)
        {
            if (_plotView.ActualModel.PlotArea.Contains(targetPoint.X, targetPoint.Y))
            {
                _plotView.UpdateLayout();
                var plotCoors = Screen2PlotCoords(targetPoint);
                try
                {
                    SelectedControlPoint.ControlPoint.SetValues(plotCoors.X, plotCoors.Y);
                }
                catch (OverflowException)
                {
                    StopActions();
                    SelectedControlPoint.ControlPoint.SetValues(plotCoors.X, plotCoors.Y);
                }
            }
        }
        private void FillWithClosestPoints(DataPoint dataPoint, Dictionary<DataPoint, double> pointsArray, Line2DOxyController lineController)
        {
            var dist = lineController.GetShortestDistance(dataPoint, IsEndPointBindingOn, IsNearPointBindingOn, IsControlPointBindingOn && lineController.IsSelected);
            if (dist != null && dist.Item1 < Double.PositiveInfinity && !pointsArray.ContainsKey(dist.Item2))
                pointsArray.Add(dist.Item2, dist.Item1);
        }

        private DataPoint Screen2PlotCoords(ScreenPoint spnt)
        {
            return _plotView.ActualModel.Axes[0].InverseTransform(spnt.X, spnt.Y, _plotView.ActualModel.Axes[1]);
        }
        private ScreenPoint PlotCoords2ScenePoint(DataPoint dPoint)
        {
            return _plotView.ActualModel.Axes[0].Transform(dPoint.X, dPoint.Y, _plotView.ActualModel.Axes[1]);
        }

        private void InitSceneObject()
        {
            var oe = _editedObject.OpticalEntity;
            if (oe == null) return;

            IsSpeciallyEditable = oe.IsSpeciallyEditable;

            if (!InitOpticalEntity(oe))
                foreach (var entity in oe.Children)
                    InitOpticalEntity(entity);

            // Init Active Lines
            foreach (var axisymSurface in _activeSurfaces)
            {
                var cLine = axisymSurface.ProfileLine2D;
                var cController = Line2DOxyController.GetLineController(cLine, _surfaceToCSMap[axisymSurface]);

                _lineToSurfaceMap[cLine] = axisymSurface;
                _lineToControllerMap[cLine] = cController;
                InitLineController(cController);
            }
            if (_lineToControllerMap.Count > 0)
                SelectedLine = _lineToControllerMap.First().Value;

            // Init Passive Lines
            foreach (var surface in _passiveSurfaces)
            {
                var cLine = new SegLine2D(surface.GetProfilePoints2D());
                var cController = Line2DOxyController.GetLineController(cLine, _surfaceToCSMap[surface]);
                cController.IsEnabled = false;

                _passiveSurfaceToControllerMap[surface] = cController;
                InitLineController(cController);
            }

            if (IsSolid)
            {
                _activeSurfacesLinks = AxisymSurfaceBinder.CreateLinks(_lineToSurfaceMap.Keys.ToList());
                foreach (var surface in _activeSurfaces)
                {
                    var curLine = surface.ProfileLine2D;

                    foreach (var link in _activeSurfacesLinks)
                        foreach (var pair in link.LinkedControlPoints)
                            if (pair.Value == curLine) _lineToControllerMap[curLine].OverlapedControlPoints.Add(pair.Key);
                }
            }

            _plotView.UpdateLayout();
        }
        private void InitLineController(Line2DOxyController cController)
        {
            foreach (var line in cController.OxyPlotSeries)
                _plotView.Series.Add(line);
            cController.OxyPlotSeries.CollectionChanged += (s, e) =>
            {
                switch (e.Action)
                {
                    case NotifyCollectionChangedAction.Add:
                        foreach (LineSeries line in e.NewItems)
                            _plotView.Series.Add(line);
                        break;
                    case NotifyCollectionChangedAction.Remove:
                        foreach (LineSeries line in e.OldItems)
                            _plotView.Series.Remove(line);
                        break;
                }
            };

            cController.PropertyChanged += (o, e) =>
                {
                    if (e.PropertyName.Equals("ControlPoints") && SelectedLine == cController)
                        SelectedControlPoint = null;
                };
        }
        private bool InitOpticalEntity(OpticalEntity oe)
        {
            if (oe.IsSurface)
            {
                AxisymLine2DSurface surf;
                if (IsSpeciallyEditable) surf = oe.OEOrigin as AxisymBezierSurface;
                else surf = oe.OEOrigin as AxisymLine2DSurface;
                if (surf == null)
                    _passiveSurfaces.Add((Surface)oe.OEOrigin);
                else
                    _activeSurfaces.Add(surf);

                _surfaceToCSMap.Add((Surface)oe.OEOrigin, oe.CS);

                return true;
            }

            return false;
        }

        #region Events

        [field: NonSerialized] public event EventHandler<ObjectChangedEventArgs> ObjectChanged;
        [field: NonSerialized] public event PropertyChangedEventHandler PropertyChanged;

        protected virtual void OnPropertyChanged([CallerMemberName] string propertyName = null)
        {
            var handler = PropertyChanged;
            if (handler != null) handler(this, new PropertyChangedEventArgs(propertyName));
        }
        protected void OnObjectChanged(EventIDs id)
        {
            var ocargs = new ObjectChangedEventArgs() { Source = this, ID = id };

            if (ObjectChanged != null)
                ObjectChanged(this, ocargs);
        }

        public enum EventIDs { UpdatePlot }

        #endregion Events
    }
}
