using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Runtime.CompilerServices;
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
    internal struct Line2DOxyHitTestResult
    {
        public bool IsLineHitted;
        public bool IsControlPointHitted;
        public Line2DOxyController HittedLine;
        public ControlPoint2DOxyController HittedControlPoint;
        public ScreenPoint Target;
        public ScreenPoint ClosestLinePoint;
        public ScreenPoint ClosestControlPoint;

        public Line2DOxyHitTestResult(ScreenPoint target, OxyPlot.HitTestResult lineHitResult, OxyPlot.HitTestResult pointHitResult)
        {
            IsLineHitted = false;
            IsControlPointHitted = false;
            Target = target;
            ClosestLinePoint = ScreenPoint.Undefined;
            ClosestControlPoint = ScreenPoint.Undefined;
            HittedControlPoint = null;
            HittedLine = null;

            if (lineHitResult != null)
            {
                IsLineHitted = true;
                ClosestLinePoint = lineHitResult.NearestHitPoint;
            }

            if (pointHitResult != null)
            {
                IsControlPointHitted = true;
                ClosestControlPoint = pointHitResult.NearestHitPoint;
            }
        }

        public double Distance2Line()
        {
            if (!IsLineHitted) return Double.PositiveInfinity;
            return (Target - ClosestLinePoint).Length;
        }

        public double Distance2ControlPoint()
        {
            if (!IsControlPointHitted) return Double.PositiveInfinity;
            return (Target - ClosestControlPoint).Length;
        }
    }

    internal abstract class Line2DOxyController : INotifyPropertyChanged
    {
        private ProfileEditorMode _displayMode;
        private bool _isEnabled;
        private bool _isSelected;
        private LineSeries _oxyLine;
        protected Dictionary<ControlPoint2D, ControlPoint2DOxyController> _point2ControllerMap;
        private SegLine2D _linePoints;
        private SegLine2D _lineControlPoints;

        public ObservableCollection<LineSeries> OxyPlotSeries { get; private set; }

        public bool IsSelected
        {
            get { return _isSelected; }
            set
            {
                if (value == IsSelected) return;
                _isSelected = value;
                UpdateLayout();
                OnPropertyChanged();
            }
        }

        public bool IsEnabled{ get { return _isEnabled; }
            set
            {
                if (!value) IsSelected = false;
                _isEnabled = value;
                UpdateLayout();
                OnPropertyChanged();
            }
        }

        internal abstract Line2D Line2D { get; }

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
                UpdateGeometry();
                UpdateLayout();
                foreach (var pair in _point2ControllerMap)
                    pair.Value.DisplayMode = value;
            }
        }

        internal CoordinateSystem CastCS { get; private set; }

        internal ObservableCollection<ControlPoint2D> OverlapedControlPoints { get; private set; }

        internal Line2DOxyHitTestResult HitLineTest(HitTestArguments args)
        {
            var lineHitResult = _oxyLine.InternalSeries.HitTest(args);
            var pointsHitResults = new SortedDictionary<double, Tuple<ControlPoint2DOxyController, OxyPlot.HitTestResult>>();
            foreach (var element in _point2ControllerMap)
            {
                var curResult = element.Value.PointSeries.InternalSeries.HitTest(args);
                if (curResult != null)
                    pointsHitResults.Add((curResult.NearestHitPoint - args.Point).Length, new Tuple<ControlPoint2DOxyController, OxyPlot.HitTestResult>(element.Value, curResult));
            }

            return new Line2DOxyHitTestResult(args.Point, lineHitResult, (pointsHitResults.Count > 0) ? pointsHitResults.First().Value.Item2 : null)
            {
                HittedControlPoint = (pointsHitResults.Count > 0) ? pointsHitResults.First().Value.Item1 : null,
                HittedLine = this
            };
        }

        internal void RemoveControlPoint(ControlPoint2DOxyController closestPoint)
        {
            Line2D.RemoveControlPoint(closestPoint.ControlPoint);
        }

        internal void AddControlPoint(DataPoint dataPoint)
        {
            Line2D.AddControlPoint(GetLocalPoint2D(dataPoint));
        }

        internal Tuple<double, DataPoint> GetShortestDistance(DataPoint curPnt, bool endPoints, bool nearPoints, bool controlPoints)
        {
            var localPoint = GetLocalPoint2D(curPnt);
            if (LinePoints.NumPoints <= 0) return null;

            Distance lineDistance = Distance.InfPoint, cPointsDistance = Distance.InfPoint;

            if (endPoints || nearPoints)
                lineDistance = LinePoints.GetShortestDistance(localPoint, endPoints, nearPoints, false);
            if (controlPoints)
                cPointsDistance = LineControlPoints.GetShortestDistance(localPoint, false, false, true);

            var resDist = Distance.GetClosest(lineDistance, cPointsDistance);
            return new Tuple<double, DataPoint>(resDist.Value, GetGlobalDataPoint(resDist.Point));
        }

        protected Line2DOxyController(CoordinateSystem castCS)
        {
            OverlapedControlPoints = new ObservableCollection<ControlPoint2D>();
            _point2ControllerMap = new Dictionary<ControlPoint2D, ControlPoint2DOxyController>();

            _oxyLine = new LineSeries();
            OxyPlotSeries = new ObservableCollection<LineSeries>() { _oxyLine };
            
            CastCS = castCS;
            _isSelected = false;
            _isEnabled = true;
            _displayMode = ProfileEditorMode.Select;
        }

        protected SegLine2D LinePoints
        {
            get
            {
                if (_linePoints == null && Line2D != null)
                {
                    _linePoints = new SegLine2D(Line2D.GetLinePoints2D());
                }
                return _linePoints;
            }
        }
        protected SegLine2D LineControlPoints
        {
            get
            {
                if (_lineControlPoints == null && Line2D != null)
                {
                    var pnts = new Points2D(Line2D.ControlPoints.Count);
                    for (int i = 0; i < Line2D.ControlPoints.Count; i++)
                        pnts[i] = Line2D.ControlPoints[i].Point;
                    _lineControlPoints = new SegLine2D(pnts);
                }
                return _lineControlPoints;
            }
        }

        protected virtual void UpdateGeometry()
        {
            _linePoints = null;
            _lineControlPoints = null;

            var line = new List<DataPoint>();
            for (int i = 0; i < LinePoints.NumPoints; i++)
                line.Add(GetGlobalDataPoint(LinePoints.Curve[i]));

            _oxyLine.ItemsSource = line;
        }

        protected virtual void UpdateLayout()
        {
            _oxyLine.LineStyle = LineStyle.Solid;

            if (!IsEnabled)
            {
                _oxyLine.LineStyle = LineStyle.Dash;
                _oxyLine.StrokeThickness = ResourceController.LineWidth;
                _oxyLine.Color = ResourceController.DisabledColor;
                return;
            }

            switch (DisplayMode)
            {
                case ProfileEditorMode.OptimizeNode:
                    _oxyLine.StrokeThickness = ResourceController.LineWidth;
                    _oxyLine.Color = ResourceController.DisabledColor;
                    break;
                case ProfileEditorMode.LockNode:
                    _oxyLine.StrokeThickness = ResourceController.LineWidth;
                    _oxyLine.Color = ResourceController.DisabledColor;
                    break;
                case ProfileEditorMode.Segline2Spline:
                case ProfileEditorMode.Spline2Segline:
                case ProfileEditorMode.CutLine:
                    _oxyLine.StrokeThickness = ResourceController.SelectedLineWidth;
                    _oxyLine.Color = ResourceController.SelectedLineColor;
                    break;
                case ProfileEditorMode.View:
                    _oxyLine.StrokeThickness = ResourceController.LineWidth;
                    _oxyLine.Color = ResourceController.LineColor;
                    break;
                default:
                    if (IsSelected)
                    {
                        _oxyLine.StrokeThickness = ResourceController.SelectedLineWidth;
                        _oxyLine.Color = ResourceController.SelectedLineColor;
                    }
                    else
                    {
                        _oxyLine.StrokeThickness = ResourceController.LineWidth;
                        _oxyLine.Color = ResourceController.LineColor;
                    }
                    break;
            }
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

        protected void InitControlPoints()
        {
            foreach (var controlPoint2D in Line2D.ControlPoints)
            {
                _point2ControllerMap[controlPoint2D] = ControlPoint2DOxyController.GetPointController(this, controlPoint2D, CastCS);
                OxyPlotSeries.Add(_point2ControllerMap[controlPoint2D].PointSeries);
            }
        }

        protected void RefreshControlPoints()
        {
            foreach (var pair in _point2ControllerMap)
            {
                OxyPlotSeries.Remove(pair.Value.PointSeries);
            }
            _point2ControllerMap.Clear();

            InitControlPoints();
        }

        internal static Line2DOxyController GetLineController(Line2D line, CoordinateSystem cs)
        {
            var bLine = line as BezierLine2D;
            if (bLine != null)
                return new BezierLine2DOxyController(bLine, cs);

            var sLine = line as SegLine2D;
            if (sLine != null)
                return new SegLine2DOxyController(sLine, cs);

            throw new ArgumentException();
        }

        #region Events

        [field: NonSerialized]
        public event PropertyChangedEventHandler PropertyChanged;

        protected virtual void OnPropertyChanged([CallerMemberName] string propertyName = null)
        {
            var handler = PropertyChanged;
            if (handler != null) handler(this, new PropertyChangedEventArgs(propertyName));
        }
        
        #endregion Events
    }
}
