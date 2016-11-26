using System;
using System.Numerics;

namespace Voronoi
{
    using TFloat = Single;
    using TNumber = Single;
    using TVector2 = Vector2;

    public static class BeachlineNodeOperations
    {
        public static Node FindDeepestNodeAtX(this Node root, TNumber x, TNumber sweepY)
        {
            var current = root;
            while (!current.IsLeaf())
            {
                var ray = (EdgeRayData)current.Data;
                var parabola = (ParabolaData)current.Left.GetRightmost().Data;
                var point = MathUtil.ComputeRayParabolaIntersection(ray.Origin, ray.Direction, parabola.Focus, sweepY);
                current = x < point.X ? current.Left : current.Right;
            }
            return current;
        }

        public static Node ComputeThreeParabolasFromDifferentYParabolaNodeCut(
            TVector2 cutParabolaFocus,
            TVector2 newParabolaFocus,
            out Node leftParabola,
            out Node centerParabola,
            out Node rightParabola
        ) {
            var intersect = MathUtil.ComputeParabolaPointGivenX(newParabolaFocus.X, newParabolaFocus.Y, cutParabolaFocus);
            var cutParabolaFocusToNewParabolaFocus = newParabolaFocus - cutParabolaFocus; // this will have a positive y, no guarantee on x sign.
            var leftEdgeDirection = new TVector2(-cutParabolaFocusToNewParabolaFocus.Y, cutParabolaFocusToNewParabolaFocus.X);
            var rightEdgeDirection = new TVector2(cutParabolaFocusToNewParabolaFocus.Y, -cutParabolaFocusToNewParabolaFocus.X);

            return new Node(new EdgeRayData(true, intersect, leftEdgeDirection))
            {
                Left = leftParabola = new Node(new ParabolaData(cutParabolaFocus)),
                Right = new Node(new EdgeRayData(true, intersect, rightEdgeDirection))
                {
                    Left = centerParabola = new Node(new ParabolaData(newParabolaFocus)),
                    Right = rightParabola = new Node(new ParabolaData(cutParabolaFocus))
                }
            };
        }
    }
}