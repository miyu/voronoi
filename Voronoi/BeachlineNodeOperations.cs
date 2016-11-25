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
                var data = (EdgeRayData)current.Data;
                var point = MathUtil.ComputeRayPointAtY(data.Origin, data.Direction, sweepY);
                if (point.X == x)
                {
                    return current;
                }
                current = x < point.X ? current.Left : current.Right;
            }
            return current;
        }

        public static Node ComputeSameYParabolaNodeCut(TVector2 cutParabolaFocus, TVector2 newParabolaFocus, out Node leftParabola, out Node rightParabola)
        {
            var leftX = cutParabolaFocus.X < newParabolaFocus.X ? cutParabolaFocus.X : newParabolaFocus.X;
            var rightX = cutParabolaFocus.X < newParabolaFocus.X ? newParabolaFocus.X : cutParabolaFocus.X;

            var middleX = leftX + (rightX - leftX) / 2;
            var between = new TVector2(middleX, newParabolaFocus.Y);
            var down = new TVector2(0, 1);
            return new Node(new EdgeRayData(between, down))
            {
                Left = leftParabola = new Node(new ParabolaData(new TVector2(leftX, newParabolaFocus.Y))),
                Right = rightParabola = new Node(new ParabolaData(new TVector2(rightX, newParabolaFocus.Y)))
            };
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

            return new Node(new EdgeRayData(intersect, leftEdgeDirection))
            {
                Left = leftParabola = new Node(new ParabolaData(cutParabolaFocus)),
                Right = new Node(new EdgeRayData(intersect, rightEdgeDirection))
                {
                    Left = centerParabola = new Node(new ParabolaData(newParabolaFocus)),
                    Right = rightParabola = new Node(new ParabolaData(cutParabolaFocus))
                }
            };
        }
    }
}