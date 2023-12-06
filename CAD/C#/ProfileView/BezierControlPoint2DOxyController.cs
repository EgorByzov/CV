using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using MathKernel;
using ObjectModel.Optimization.ControlPoints;

namespace ViewModelCore.OpticalSceneEntities.ProfileView
{
    internal class BezierControlPoint2DOxyController : ControlPoint2DOxyController
    {
        private BezierControlPoint2D BezierCP { get; set; }

        public BezierControlPoint2DOxyController(BezierControlPoint2D cPoint, Line2DOxyController parentLine, CoordinateSystem castCS)
        {
            InitConstructor(cPoint, parentLine, castCS);
            BezierCP = cPoint;
            UpdateLayout();
        }

        protected override void UpdateLayout()
        {
            base.UpdateLayout();

            switch (DisplayMode)
            {
                case ProfileEditorMode.OptimizeNode:
                case ProfileEditorMode.LockNode:
                    break;
                case ProfileEditorMode.AddNode:
                case ProfileEditorMode.RemoveNode:
                case ProfileEditorMode.CutLine:
                    if (BezierCP.IsBaseControlPoint)
                        PointSeries.MarkerType = OxyPlot.MarkerType.Square;
                    else
                        HidePoint();
                    break;
                default:
                    if (BezierCP.IsBaseControlPoint)
                        PointSeries.MarkerType = OxyPlot.MarkerType.Square;
                    else
                        PointSeries.MarkerType = OxyPlot.MarkerType.Circle;
                    break;
            }
        }

    }
}
