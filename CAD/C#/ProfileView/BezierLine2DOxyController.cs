using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Common;
using MathKernel;
using MathKernel.VectorAlgebra;
using ObjectModel.OpticalSceneEntities.OriginEntities.Lines;
using ObjectModel.Optimization.ControlPoints;
using OxyPlot;
using OxyPlot.Wpf;
using ViewModelCore.Controllers;

namespace ViewModelCore.OpticalSceneEntities.ProfileView
{
    internal sealed class BezierLine2DOxyController : Line2DOxyController
    {
        private LineSeries _oxyControlPointContour;
        private List<DataPoint> _supportContour;
 
        public BezierLine2D BezierLine { get; private set; }

        internal override Line2D Line2D { get { return BezierLine; } }

        private bool IsContourHidden { get { return _oxyControlPointContour.ItemsSource == null; } }

        internal BezierLine2DOxyController(BezierLine2D line, CoordinateSystem castCS) : base(castCS)
        {
            BezierLine = line;

            _supportContour = new List<DataPoint>();
            _oxyControlPointContour = new LineSeries();

            OxyPlotSeries.Add(_oxyControlPointContour);

            System.Windows.WeakEventManager<BezierLine2D, ObjectChangedEventArgs>.AddHandler(BezierLine, "ObjectChanged",
                (o, e) =>
                {
                    if (e.ID.Equals(Line2D.EventIDs.Geometry))
                        UpdateGeometry();
                    else if (e.ID.Equals(Line2D.EventIDs.ControlPoints))
                    {
                        OnPropertyChanged("ControlPoints");
                        RefreshControlPoints();
                        UpdateGeometry();
                        UpdateLayout();
                    }
                });

            InitControlPoints();
            UpdateGeometry();
            UpdateLayout();
        }

        private void HideContour()
        {
            _oxyControlPointContour.ItemsSource = null;
        }

        private void ShowContour()
        {
            _oxyControlPointContour.ItemsSource = _supportContour;
        }

        protected override void UpdateGeometry()
        {
            base.UpdateGeometry();

            _supportContour.Clear();
            var CPoints = BezierLine.BezierControlPoints;
            foreach (var curPoint in CPoints)
                _supportContour.Add(GetGlobalDataPoint(curPoint.Point));

            if (!IsContourHidden) ShowContour();
        }

        protected override void UpdateLayout()
        {
            base.UpdateLayout();

            HideContour();
            if (!IsEnabled) return;

            switch (DisplayMode)
            {
                case ProfileEditorMode.AddNode:
                case ProfileEditorMode.RemoveNode:
                case ProfileEditorMode.CutLine:
                    break;
                case ProfileEditorMode.OptimizeNode:
                case ProfileEditorMode.LockNode:
                case ProfileEditorMode.Segline2Spline:
                case ProfileEditorMode.Spline2Segline:                
                case ProfileEditorMode.View:
                    _oxyControlPointContour.StrokeThickness = ResourceController.ContourWidth;
                    _oxyControlPointContour.Color = ResourceController.DisabledColor;
                    ShowContour();
                    break;
                default:
                    if (IsSelected)
                    {
                        ShowContour();
                        _oxyControlPointContour.StrokeThickness = ResourceController.ContourWidth;
                        _oxyControlPointContour.Color = ResourceController.ContourLineColor;
                    }
                    break;
            }
        }
    }
}
