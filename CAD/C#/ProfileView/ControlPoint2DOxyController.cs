using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;
using MathKernel;
using MathKernel.VectorAlgebra;
using ObjectModel.OpticalSceneEntities.OriginEntities.Lines;
using ObjectModel.Optimization.ControlPoints;
using OxyPlot;
using OxyPlot.Wpf;
using ViewModelCore.Controllers;

namespace ViewModelCore.OpticalSceneEntities.ProfileView
{
    internal class ControlPoint2DOxyController
    {
        private ProfileEditorMode _displayMode;
        private LineSeries _oxySeries;
        private bool _isSelected;
        private DataPoint _supportPoint;

        internal Line2DOxyController ParentLineController { get; private set; }
        internal ControlPoint2D ControlPoint { get; private set; }
        internal DataPoint OxyPoint {get { return _supportPoint; }}
        internal CoordinateSystem CastCS { get; private set; }
        internal bool IsOverlaped { get; private set; }
        internal ProfileEditorMode DisplayMode
        {
            get
            {
                return _displayMode;
            }
            set
            {
                if (DisplayMode == value) return;

                _displayMode = value;
                UpdateLayout();
            }
        }

        public LineSeries PointSeries { get { return _oxySeries; } private set { _oxySeries = value; } }
        public bool IsSelected
        {
            get { return _isSelected; }
            set
            {
                if (IsSelected == value) return;
                _isSelected = value;
                UpdateLayout();
            }
        }
        public bool IsLocked { get { return ControlPoint.IsLocked; } }
        public bool IsOptimized { get { return ControlPoint.IsOptimizationOn; } }
        public bool IsUnlockAvailable { get { return ControlPoint.IsUnlockAllowed; } }
        protected bool IsHidden { get { return PointSeries.ItemsSource == null; } }

        protected ControlPoint2DOxyController(){}
        public ControlPoint2DOxyController(ControlPoint2D cPoint, Line2DOxyController parentLine, CoordinateSystem castCS)
        {
            InitConstructor(cPoint, parentLine, castCS);
            UpdateLayout();
        }
        
        internal double HitPointTest(HitTestArguments args)
        {
            return Double.PositiveInfinity;
        }

        protected virtual void InitConstructor(ControlPoint2D cPoint, Line2DOxyController parentLine, CoordinateSystem castCS)
        {
            ControlPoint = cPoint;
            ParentLineController = parentLine;
            CastCS = castCS;

            PointSeries = new LineSeries();

            // Bind Control Point
            cPoint.ObjectChanged += (sender, args) =>
            {
                if (args.ID.Equals(ObjectModel.Optimization.ControlPoints.ControlPoint.EventIDs.Geometry) || 
                    args.ID.Equals(ObjectModel.Optimization.ControlPoints.ControlPoint.EventIDs.PointObject))
                    UpdateGeometry();
                else
                    UpdateLayout();
            };

            // Bind Parent Line
            ParentLineController.OverlapedControlPoints.CollectionChanged += (sender, args) =>
            {
                if (args.Action == NotifyCollectionChangedAction.Add)
                    foreach (ControlPoint2D newPoint in args.NewItems)
                        if (newPoint == ControlPoint) IsOverlaped = true;
                        else if (args.Action == NotifyCollectionChangedAction.Remove)
                            foreach (ControlPoint2D oldPoint in args.OldItems)
                                if (oldPoint == ControlPoint) IsOverlaped = false;
            };
            foreach (ControlPoint2D controlPoint in ParentLineController.OverlapedControlPoints)
                if (controlPoint == ControlPoint) IsOverlaped = true;

            System.Windows.WeakEventManager<Line2DOxyController, PropertyChangedEventArgs>.AddHandler(parentLine, "PropertyChanged",
                (o, e) =>
                {
                    UpdateLayout();
                });

            UpdateGeometry();
        }

        protected virtual void UpdateGeometry()
        {
            _supportPoint = GetGlobalDataPoint(ControlPoint.Point);
            if (!IsHidden) ShowPoint();
        }

        protected virtual void UpdateLayout()
        {
            PointSeries.LineStyle = LineStyle.None;
            PointSeries.CanTrackerInterpolatePoints = false;            

            if (!ParentLineController.IsEnabled || DisplayMode == ProfileEditorMode.View)
            {
                HidePoint();
                return;
            }
            ShowPoint();

            switch (DisplayMode)
            {
                case ProfileEditorMode.OptimizeNode:
                    PointSeries.MarkerType = MarkerType.Circle;
                    PointSeries.MarkerSize = ResourceController.SelectedPointSize;
                    PointSeries.MarkerStroke = ResourceController.PointColor;

                    if (ControlPoint.IsLocked)
                        PointSeries.MarkerFill = ResourceController.DisabledColor;
                    else if (IsOptimized)
                        PointSeries.MarkerFill = ResourceController.OptimizedPointColor;
                    else
                        PointSeries.MarkerFill = Colors.Transparent;

                    break;
                case ProfileEditorMode.LockNode:
                    PointSeries.MarkerType = MarkerType.Circle;
                    PointSeries.MarkerSize = ResourceController.SelectedPointSize;
                    PointSeries.MarkerStroke = ResourceController.PointColor;

                    if (ControlPoint.IsLocked)
                        PointSeries.MarkerFill = ResourceController.DisabledColor;
                    else
                        PointSeries.MarkerFill = Colors.Transparent;

                    break;
                case ProfileEditorMode.Segline2Spline:
                case ProfileEditorMode.Spline2Segline:
                case ProfileEditorMode.CutLine:
                    PointSeries.MarkerType = MarkerType.Circle;
                    PointSeries.MarkerSize = ResourceController.PointSize;
                    PointSeries.MarkerFill = ResourceController.PointColor;
                    PointSeries.MarkerStroke = ResourceController.PointColor;
                    break;
                default:
                    if (!ParentLineController.IsSelected)
                    {
                        HidePoint();
                        return;
                    }
                    PointSeries.MarkerStroke = Colors.Transparent;
                    PointSeries.MarkerType = MarkerType.Square;
                    if (IsSelected)
                    {
                        PointSeries.MarkerSize = ResourceController.SelectedPointSize;
                        PointSeries.MarkerFill = ResourceController.SelectedPointColor;
                        PointSeries.MarkerStroke = ResourceController.ContourLineColor;
                    }
                    else
                    {
                        PointSeries.MarkerSize = ResourceController.PointSize;
                        PointSeries.MarkerFill = ResourceController.PointColor;
                    }
                    break;
            }
        }

        protected void HidePoint()
        {
            PointSeries.ItemsSource = null;
        }

        protected void ShowPoint()
        {
            PointSeries.ItemsSource = new List<DataPoint>() { _supportPoint };
        }

        protected DataPoint GetGlobalDataPoint(Point2D lPoint)
        {
            var gPoint = CastCS.GetGlobalPoint(new Point(lPoint.X, 0, lPoint.Y));
            return new DataPoint(gPoint.X, gPoint.Z);
        }

        protected Point2D GetLocalPoint2D(DataPoint gPoint)
        {
            var lPoint = CastCS.GetLocalPoint(new Point(gPoint.X, 0, gPoint.Y));
            return new Point2D(lPoint.X, lPoint.Z);
        }

        internal static ControlPoint2DOxyController GetPointController(Line2DOxyController parentLine, ControlPoint2D point, CoordinateSystem cs)
        {
            var bPoint = point as BezierControlPoint2D;
            if (bPoint != null)
                return new BezierControlPoint2DOxyController(bPoint, parentLine, cs) { DisplayMode = parentLine.DisplayMode };

            return new ControlPoint2DOxyController(point, parentLine, cs) {DisplayMode = parentLine.DisplayMode};
        }

    }
}
