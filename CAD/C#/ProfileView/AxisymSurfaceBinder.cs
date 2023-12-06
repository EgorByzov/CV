using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathKernel.VectorAlgebra;
using ObjectModel.Functions.Fitting;
using ObjectModel.OpticalSceneEntities;
using ObjectModel.OpticalSceneEntities.OriginEntities.Lines;
using ObjectModel.OpticalSceneEntities.OriginEntities.Surfaces;
using ObjectModel.Optimization.ControlPoints;

namespace ViewModelCore.OpticalSceneEntities.ProfileView
{
    public class AxisymSurfaceBinder
    {
        internal ControlPoint2D GeneralizedControlPoint { get; private set; }

        internal Line2D GeneralizedLine2D { get; private set; }

        internal Dictionary<ControlPoint2D, Line2D> LinkedControlPoints { get; private set; }

        private AxisymSurfaceBinder(ControlPoint2D basePoint, Line2D baseLine, Dictionary<ControlPoint2D, Line2D> linkedLines)
        {
            GeneralizedControlPoint = basePoint;
            GeneralizedLine2D = baseLine;
            LinkedControlPoints = linkedLines;
            InitLinks();
        }

        private void InitLinks()
        {
        }

        internal static List<AxisymSurfaceBinder> CreateLinks(List<Line2D> lines)
        {
            var result = new List<AxisymSurfaceBinder>();

            var pointColumn = new Dictionary<ControlPoint2D, Line2D>();

            for (int i = 0; i < lines.Count; i++)
            {
                var curLine = lines[i];
                pointColumn[curLine.FirstControlPoint] = curLine;
                pointColumn[curLine.LastControlPoint] = curLine;
            }


            while (pointColumn.Count > 0)
            {
                var curPoint = pointColumn.First().Key;
                var curLine = pointColumn.First().Value;
                var curLinkedPoints = new Dictionary<ControlPoint2D, Line2D>();
                pointColumn.Remove(curPoint);
                
                foreach (var keyValuePair in pointColumn)
                {
                    var linkedPoint = keyValuePair.Key;
                    var linkedLine = keyValuePair.Value;
                    if (curPoint.Point.Equals(linkedPoint.Point))
                    {
                        curLinkedPoints[linkedPoint] = linkedLine;
                    }
                }

                if (curLinkedPoints.Count > 0)
                {
                    foreach (var item in curLinkedPoints)
                        pointColumn.Remove(item.Key);

                    result.Add(new AxisymSurfaceBinder(curPoint, curLine, curLinkedPoints));
                }
            }

            return result;
        }
    }
}
