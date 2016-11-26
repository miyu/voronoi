using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Threading;
using System.Windows.Forms;

namespace Voronoi
{
//    using TNumber = Int64;
//    using TVector2 = IntVector2;
    using TFloat = Single;
    using TNumber = Single;
    using TVector2 = Vector2;

    public static class MathUtil
    {
        public static TNumber Average(TNumber a, TNumber b)
        {
            if (a > b)
            {
                return Average(b, a);
            }
            return a + (b - a) / 2;
        }

        public static bool ApproximatelyEqual(TFloat a, TFloat b)
        {
            const TFloat epsilon = 1E-4f;
            var diff = Math.Abs(a - b);
            // ReSharper disable CompareOfFloatsByEqualityOperator
            return a == b
                   || (a == 0 || b == 0 || diff < TFloat.Epsilon
                       ? diff < TFloat.Epsilon
                       : diff / (Math.Abs(a) + Math.Abs(b)) < epsilon);
            // ReSharper restore CompareOfFloatsByEqualityOperator
        }

        public static TVector2 ComputeParabolaPointGivenX(TNumber rx, TNumber directrixY, TVector2 focus)
        {
            // Let result be a point r (rx, ry) located between the focus's y and directrix y.
            // directrixY - ry = sqrt((rx - focusX)^2 + (ry - focusY)^2)
            // (directrixY - ry)^2 = (rx - focusX)^2 + (ry - focusY)^2
            // directrixY^2 - 2 directrixY ry + ry^2 = rx^2 - 2 rx focusX + focusX^2 + ry^2 - 2 ry focusY + focusY^2
            // -2 directrixY ry + ry^2 = rx^2 - 2 rx focusX + focusX^2 + ry^2 - 2 ry focusY + focusY^2 - directrixY^2
            // -2 directrixY ry + ry^2 - ry^2 + 2 ry focusY = rx^2 - 2 rx focusX + focusX^2 + focusY^2 - directrixY^2
            // -2 directrixY ry + 2 ry focusY = rx^2 - 2 rx focusX + focusX^2 + focusY^2 - directrixY^2
            // ry (-2 directrixY + 2 focusY) = rx^2 - 2 rx focusX + focusX^2 + focusY^2 - directrixY^2
            // ry = (rx^2 - 2 rx focusX + focusX^2 + focusY^2 - directrixY^2) / (-2 directrixY + 2 focusY)
            var numerator = rx * rx - 2 * rx * focus.X + focus.X * focus.X + focus.Y * focus.Y - directrixY * directrixY;
            var denominator = -2 * directrixY + 2 * focus.Y;
            var ry = numerator / denominator;
            return new TVector2(rx, ry);
        }

        private static void ComputeParabolaStandardForm(
            TVector2 focus, TNumber directrixY,
            out TNumber a, out TNumber b, out TNumber c)
        {
            if (ApproximatelyEqual(focus.Y, directrixY))
            {
                a = b = 0;
                c = directrixY;
                return;
            }
            var denominator = -2 * directrixY + 2 * focus.Y;
            a = 1.0f / denominator;
            b = -2 * focus.X / denominator;
            c = (focus.X * focus.X + focus.Y * focus.Y - directrixY * directrixY) / denominator;
        }

        public static void ComputeParabolaIntersections(
            TVector2 focusA, TVector2 focusB, TNumber directrixY,
            out TVector2 leftIntersection, out TVector2 rightIntersection)
        {
            // Let a and b the foci of two parabolae (the near-pretentious plural of parabola).
            // Let result be a point r (rx, ry) located at equal distances between the foci and directrix.
            // Of course, with two parabolas it's easy to see there can be two such intersections.
            //
            // First convert to standard form (y = a x^2 + bx + c).
            //    Let p be a point (px, py) on the parabola.
            //    Let f be a point (fx, fy) representing the parabola's focus
            //    Let dy be the y-coordinate of the directrix.
            //
            //        dy - py = |p - f|
            //        dy - py = sqrt((px - fx)^2 + (py - fy)^2)
            //        (dy - py)^2 = (px - fx)^2 + (py - fy)^2
            //        dy^2 - 2 dy py + py^2 = px^2 - 2 px fx + fx^2 + py^2 - 2 py fy + fy^2
            //        dy^2 - 2 dy py = px^2 - 2 px fx + fx^2 - 2 py fy + fy^2
            //        - 2 dy py + 2 py fy = px^2 - 2 px fx + fx^2 + fy^2 - dy^2
            //        py (-2dy + 2fy) = px^2 - 2 px fx + fx^2 + fy^2 - dy^2
            //        py = (px^2 - 2 px fx + fx^2 + fy^2 - dy^2) / (-2dy + 2fy)
            //        py = (1 / (-2dy + 2fy))(px^2) + (-2fx / (-2dy + 2fy)) px + (fx^2 + fy^2 - dy^2) / (-2dy + 2fy)
            //             a                  px^2  + b                     px + c
            //
            // Given two parabolas in standard form (y1 = a1 x^2 + b1 x + c1, y2 = a2 x^2 + b2 x + c2),
            // solve the system of equations...
            //
            //                        y1 = y2
            //        a1 x^2 + b1 x + c1 = a2 x^2 + b2 x + c2
            //                         0 = (a2 - a1) x^2 + (b2 - b1) x + (c2 - c1)
            //
            //                         x = -(b2 - b1) +/- sqrt((b2 - b1)^2 - 4(a2 - a1)(c2 - c1))
            //                             ------------------------------------------------------
            //                              2(a2 - a1)
            TNumber a1, b1, c1, a2, b2, c2;
            ComputeParabolaStandardForm(focusA, directrixY, out a1, out b1, out c1);
            ComputeParabolaStandardForm(focusB, directrixY, out a2, out b2, out c2);

            TNumber a = a2 - a1, b = b2 - b1, c = c2 - c1;
            var discriminant = b * b - 4 * a * c;
            var xCenter = -b / (2 * a);
            var xOffset = (TNumber)Math.Abs(Math.Sqrt(discriminant) / (2 * a));
            var leftX = xCenter - xOffset;
            var rightX = xCenter + xOffset;
            leftIntersection = new Vector2(leftX, a1 * leftX * leftX + b1 * leftX + c1);
            rightIntersection = new Vector2(rightX, a1 * rightX * rightX + b1 * rightX + c1);
        }

        public static TVector2 ComputeRayParabolaIntersection(
            TVector2 rayOrigin, TVector2 rayDirection,
            TVector2 parabolaFocus, TNumber directrixY)
        {
            if (ApproximatelyEqual(0, rayDirection.X))
            {
                var x = rayOrigin.X;
                TNumber sa, sb, sc;
                ComputeParabolaStandardForm(parabolaFocus, directrixY, out sa, out sb, out sc);
                return new Vector2(x, sa * x * x + sb * x + sc);
            }

            if (MathUtil.ApproximatelyEqual(directrixY, parabolaFocus.Y))
            {
                directrixY += 0.2f;
            }

            // Given an intersecting ray R and parabola P, find the first point of intersection.
            //
            // Let ro be the ray origin point (rox, roy)
            //     rd be the ray direction vector (rdx, rdy)
            //     pf be the parabola focus point (pfx, pfy)
            //     dy be the directrix y-coordinate
            //
            // Let t be the ray parameterization variable in [0, inf)
            //     s be the solution point (sx, sy), s = ro + rd * t
            //
            // By definition of parabola:
            //     dy - sy = |s - f|
            //     (dy - sy)^2 = |s - f|^2
            //
            // Expand (dy - sy)^2:
            //     (dy - sy)^2
            //     (dy - (roy + rdy * t))^2
            //     ((dy - roy) - rdy * t)^2
            //     (dy - roy)^2 - 2 (dy - roy) (rdy * t) + (rdy * t)^2
            //     (dy - roy)^2 - 2 (dy - roy) (rdy) t + (rdy^2) t^2
            //
            // Expand |s - f|^2:
            //     (sx - fx)^2 + (sy - fy)^2
            //     ((rox + rdx * t) - fx)^2 + ((roy + rdy * t) - fy)^2
            //
            //     Expand ((rox + rdx * t) - fx)^2:
            //         ((rdx * t) + (rox - fx))^2
            //         (rdx * t)^2 + 2 (rdx * t) (rox - fx) + (rox - fx)^2
            //         (rdx^2) t^2 + t (2 rdx (rox - fx)) + (rox - fx)^2
            //
            //     Expand ((roy + rdy * t) - fy)^2: (likewise)
            //         (rdy^2) t^2 + t (2 rdy (roy - fy)) + (roy - fy)^2
            //
            //     (rdx^2) t^2 + t (2 rdx (rox - fx)) + (rox - fx)^2 + (rdy^2) t^2 + t (2 rdy (roy - fy)) + (roy - fy)^2
            //     (rdx^2 + rdy^2) t^2 + t (2 rdx (rox - fx) + 2 rdy (roy - fy)) + (rox - fx)^2 + (roy - fy)^2
            //
            // Back to: dy - sy = |s - f|
            //     (dy - sy)^2 = |s - f|^2
            //
            //     (dy - roy)^2 - 2 (dy - roy) (rdy) t + (rdy^2) t^2 =
            //         (rdx^2 + rdy^2) t^2 + t (2 rdx (rox - fx) + 2 rdy (roy - fy)) + (rox - fx)^2 + (roy - fy)^2
            //
            //     0 = t^2 (rdx^2 + rdy^2 - rdy^2) +
            //         t   (2 rdx (rox - fx) + 2 rdy (roy - fy) + 2 (dy - roy) (rdy)) +
            //         1   ((rox - fx)^2 + (roy - fy)^2 - (dy - roy)^2)
            //
            //     0 = t^2 (rdx^2) +
            //         t   (2 rdx rox - 2 rdx fx + 2 rdy roy - 2 rdy fy + 2 rdy dy - 2 rdy roy) +
            //         1   (rox^2 - 2 rox fx + fx^2 + roy^2 - 2 roy fy + fy^2 - dy^2 + 2 dy roy - roy^2)
            //
            //     0 = t^2 (rdx^2) +
            //         t   (2 rdx rox - 2 rdx fx - 2 rdy fy + 2 rdy dy) +
            //         1   (rox^2 - 2 rox fx + fx^2 - 2 roy fy + fy^2 - dy^2 + 2 dy roy)

            var roxMinusFx = rayOrigin.X - parabolaFocus.X;
            var royMinusFy = rayOrigin.Y - parabolaFocus.Y;
            var dyMinusRoy = directrixY - rayOrigin.Y;
            var a = rayDirection.X * rayDirection.X;
            var b = 2 * (rayDirection.X * roxMinusFx
                         - rayDirection.Y * (parabolaFocus.Y - directrixY));
            var c = roxMinusFx * roxMinusFx + royMinusFy * royMinusFy - dyMinusRoy * dyMinusRoy;
            TNumber pickedT = default(TNumber);
            if (!MathUtil.ApproximatelyEqual(0, a))
            {
                var centerT = -b / (2 * a);
                var offsetT = (TNumber) Math.Sqrt(b * b - 4 * a * c) / (2 * a);
                var lowerT = centerT - offsetT;
                var upperT = centerT + offsetT;
                pickedT = lowerT >= 0 ? lowerT : upperT;
            }
            if (TNumber.IsNaN(pickedT))
            {
//                Debugger.Break();
                pickedT = 0;
            }
            return rayOrigin + pickedT * rayDirection;
        }

        public static TVector2 ComputeRayPointAtY(TVector2 rayOrigin, TVector2 rayDirection, TNumber y)
        {
            var dy = y - rayOrigin.Y;
            return rayOrigin + (rayDirection * dy) / rayDirection.Y;
        }

        public static TNumber Det3(
            TNumber m11, TNumber m12, TNumber m13,
            TNumber m21, TNumber m22, TNumber m23,
            TNumber m31, TNumber m32, TNumber m33)
        {
            return m11 * m22 * m33 + m12 * m23 * m31 + m13 * m21 * m32 - (
                       m31 * m22 * m13 + m32 * m23 * m11 + m33 * m21 * m12
                );
        }

        public static void FindCircumcircle(TVector2 p1, TVector2 p2, TVector2 p3, out TVector2 circumcenter, out TNumber radius)
        {
            var p1DotP1 = p1.X * p1.X + p1.Y * p1.Y;
            var p2DotP2 = p2.X * p2.X + p2.Y * p2.Y;
            var p3DotP3 = p3.X * p3.X + p3.Y * p3.Y;
            var a = Det3(
                p1.X, p1.Y, 1,
                p2.X, p2.Y, 1,
                p3.X, p3.Y, 1
            );
            var bx = -Det3(
                p1DotP1, p1.Y, 1,
                p2DotP2, p2.Y, 1,
                p3DotP3, p3.Y, 1
            );
            var by = Det3(
                p1DotP1, p1.X, 1,
                p2DotP2, p2.X, 1,
                p3DotP3, p3.X, 1
            );
            var c = -Det3(
                p1DotP1, p1.X, p1.Y,
                p2DotP2, p2.X, p2.Y,
                p3DotP3, p3.X, p3.Y
            );
            radius = (TNumber)Math.Sqrt(bx * bx + by * by - 4 * a * c) / (2 * Math.Abs(a));
            circumcenter = new TVector2(-bx / (2 * a), -by / (2 * a));
        }
    }

    public struct IntVector2
    {
        public IntVector2(TNumber x, TNumber y)
        {
            X = x;
            Y = y;
        }

        public TNumber X { get; }
        public TNumber Y { get; }

        public override int GetHashCode() => X.GetHashCode() ^ Y.GetHashCode();
        public override bool Equals(object obj) => obj is IntVector2 && X == ((IntVector2) obj).X && Y == ((IntVector2) obj).Y;

        public static IntVector2 operator +(IntVector2 a, IntVector2 b) => new IntVector2(a.X + b.X, a.Y + b.Y);
        public static IntVector2 operator -(IntVector2 a, IntVector2 b) => new IntVector2(a.X - b.X, a.Y - b.Y);
        public static IntVector2 operator *(IntVector2 a, TNumber b) => new IntVector2(a.X * b, a.Y * b);
        public static IntVector2 operator *(TNumber a, IntVector2 b) => new IntVector2(b.X * a, b.Y * a);
        public static IntVector2 operator /(IntVector2 a, TNumber b) => new IntVector2(a.X / b, a.Y / b);
        public static IntVector2 operator /(TNumber a, IntVector2 b) => new IntVector2(b.X / a, b.Y / a);
    }

    public class FortunesAlgorithmState
    {
        private FortunesAlgorithmState() {}

        public ISet<Vector2> InitialPoints { get; private set; }
        public Node Root;
        public PriorityQueue<ISweepEvent> EventQueue { get; set; }
        public List<VoronoiEdge> VoronoiEdges { get; } = new List<VoronoiEdge>();
        public List<VoronoiEdge> DelanayEdges { get; } = new List<VoronoiEdge>();

        public void PrintBeachlineTree(Node highlightedNode = null)
        {
            PrintBeachlineTreeHelper(Root, highlightedNode, 1);
            Console.WriteLine("");
        }

        private void PrintBeachlineTreeHelper(Node current, Node highlightedNode, int indent)
        {
            if (current == null) return;
            Console.WriteLine(new string('\t', indent) + current.Data + (highlightedNode == current ? " ***" : ""));
            PrintBeachlineTreeHelper(current.Left, highlightedNode, indent + 1);
            PrintBeachlineTreeHelper(current.Right, highlightedNode, indent + 1);
        }

        public static FortunesAlgorithmState CreateAndInitialize(ISet<TVector2> points)
        {
            var eventComparer = Comparer<ISweepEvent>.Create(CompareSweepEvents);
            var eventQueue = new PriorityQueue<ISweepEvent>(eventComparer);

            foreach (var point in points)
            {
                eventQueue.Enqueue(new SiteEvent(point));
            }

            return new FortunesAlgorithmState
            {
                InitialPoints = points,
                Root = null,
                EventQueue = eventQueue
            };
        }

        private static int CompareSweepEvents(ISweepEvent a, ISweepEvent b)
        {
            var result = a.Y.CompareTo(b.Y);
            if (result != 0) return result;
            return ((int) a.Type).CompareTo((int) b.Type);
        }
    }

    public class FortunesAlgorithmExecutor
    {
        public void ExecuteNextEvent(FortunesAlgorithmState state)
        {
            if (state.Root == null)
            {
                ProcessFirstSiteEvents(state);
                return;
            }

            ISweepEvent sweepEvent = state.EventQueue.Dequeue();
            if (sweepEvent.Type == SweepEventType.Site)
                ProcessSiteEvent(state, (SiteEvent) sweepEvent);
            else
                ProcessCircleEvent(state, (CircleEvent) sweepEvent);
        }

        public void ExecuteToCompletion(FortunesAlgorithmState state)
        {
            while (state.EventQueue.Any())
                ExecuteNextEvent(state);

            ProcessCleanup(state);
        }

        private void ProcessFirstSiteEvents(FortunesAlgorithmState state)
        {
            var initialSiteEvent = (SiteEvent) state.EventQueue.Dequeue();

            // Pull other site events with the same y-value
            var siteEvents = new List<SiteEvent> { initialSiteEvent };
            ISweepEvent sweepEvent;
            while (state.EventQueue.TryPeek(out sweepEvent) &&
                   MathUtil.ApproximatelyEqual(sweepEvent.Y, initialSiteEvent.Y))
            {
                state.EventQueue.Dequeue();
                siteEvents.Add((SiteEvent)sweepEvent);
            }

            // Order by X and divide and conquer to build beachfront tree.
            siteEvents.Sort((a, b) => a.Point.X.CompareTo(b.Point.X));
            state.Root = BuildInitialBeachFront(siteEvents, 0, siteEvents.Count);

            foreach (var pair in siteEvents.Zip(siteEvents.Skip(1), Tuple.Create))
                state.DelanayEdges.Add(new VoronoiEdge(pair.Item1.Point, true, pair.Item2.Point, true));
        }

        private Node BuildInitialBeachFront(IReadOnlyList<SiteEvent> siteEvents, int startIndexInclusive, int endIndexExclusive)
        {
            var count = endIndexExclusive - startIndexInclusive;
            Console.WriteLine(count);

            if (count == 0)
                return null;

            if (count == 1)
                return new Node(new ParabolaData(siteEvents[startIndexInclusive].Point));

            var midpoint = startIndexInclusive + count / 2;
            var left = BuildInitialBeachFront(siteEvents, startIndexInclusive, midpoint);
            var right = BuildInitialBeachFront(siteEvents, midpoint, endIndexExclusive);

            if (left == null || right == null)
                return left ?? right;

            var leftFocus = siteEvents[midpoint - 1].Point;
            var rightFocus = siteEvents[midpoint].Point;
            var midX = MathUtil.Average(leftFocus.X, rightFocus.X);
            var edgeData = new EdgeRayData(false, new TVector2(midX, leftFocus.Y), new TVector2(0, 1));
            return new Node(edgeData)
            {
                Left = left,
                Right = right
            };
        }

        public void ProcessCleanup(FortunesAlgorithmState state)
        {
            var s = new Stack<Node>();
            s.Push(state.Root);

            while (s.Any())
            {
                var current = s.Pop();
                if (current == null) continue;

                var edgeRayData = current.Data as EdgeRayData;
                if (edgeRayData == null) continue;

                s.Push(current.Left);
                s.Push(current.Right);

                state.VoronoiEdges.Add(
                    new VoronoiEdge(
                        edgeRayData.Origin, edgeRayData.IsOriginBounded,
                        edgeRayData.Origin + edgeRayData.Direction, false));

                state.Root = null;
            }
        }

        private void ProcessSiteEvent(FortunesAlgorithmState state, SiteEvent e)
        {
            var sweepY = e.Y;
            var newParabolaFocus = e.Point;
            if (state.Root == null)
            {
                state.Root = new Node(new ParabolaData(newParabolaFocus));
                return;
            }

            var cutNode = state.Root.FindDeepestNodeAtX(newParabolaFocus.X, newParabolaFocus.Y);
            if (!cutNode.IsLeaf()) throw new NotSupportedException("Cutting edge node not supported.");

            var cutParabolaFocus = ((ParabolaData)cutNode.Data).Focus;

            Node leftParabola, centerParabola, rightParabola;
            Node newNode = BeachlineNodeOperations.ComputeThreeParabolasFromDifferentYParabolaNodeCut(
                cutParabolaFocus, newParabolaFocus,
                out leftParabola, out centerParabola, out rightParabola);
            NodeOperations.ReplaceNode(cutNode, newNode, ref state.Root);

            state.DelanayEdges.Add(new VoronoiEdge(newParabolaFocus, true, cutParabolaFocus, true));

            Node leftLeftParabola;
            if (leftParabola.TryGetLeftLeaf(out leftLeftParabola))
                HandleAddCircleEvent(state, leftLeftParabola, leftParabola, centerParabola, sweepY);

            Node rightRightParabola;
            if (rightParabola.TryGetRightLeaf(out rightRightParabola))
                HandleAddCircleEvent(state, centerParabola, rightParabola, rightRightParabola, sweepY);
        }

        void HandleAddCircleEvent(
            FortunesAlgorithmState state,
            Node leftParabola, Node centerParabola, Node rightParabola,
            TNumber sweepY)
        {
            var leftParabolaFocus = ((ParabolaData) leftParabola.Data).Focus;
            var centerParabolaFocus = ((ParabolaData) centerParabola.Data).Focus;
            var rightParabolaFocus = ((ParabolaData) rightParabola.Data).Focus;

            // center gets swallowed by left/right when sweep line hits bottom of foci circumcircle
            TVector2 circumcenter;
            TNumber radius;
            MathUtil.FindCircumcircle(
                leftParabolaFocus, centerParabolaFocus, rightParabolaFocus,
                out circumcenter, out radius);

            // All three points have a circumcenter - the parabola-swallowing
            // circle event will only happen if edge rays point towards the circumcenter.
            Node leftAncestor, rightAncestor;
            centerParabola.FindDirectionalAncestors(out leftAncestor, out rightAncestor);

            var leftEdgeRayData = (EdgeRayData) leftAncestor.Data;
            var rightEdgeRayData = (EdgeRayData) rightAncestor.Data;
            var leftEdgeOriginToCircumcenter = circumcenter - leftEdgeRayData.Origin;
            var rightEdgeOriginToCircumcenter = circumcenter - rightEdgeRayData.Origin;
            if (TVector2.Dot(leftEdgeOriginToCircumcenter, leftEdgeRayData.Direction) <= 0 ||
                TVector2.Dot(rightEdgeOriginToCircumcenter, rightEdgeRayData.Direction) <= 0)
                return;

            if (circumcenter.Y + radius > sweepY)
            {
                state.EventQueue.Enqueue(
                    new CircleEvent(
                        circumcenter,
                        radius,
                        centerParabola
                    ));
            }
        }

        void ProcessCircleEvent(FortunesAlgorithmState state, CircleEvent e)
        {
            var sweepY = e.Y;
            var swallowedParabolaNode = e.SwallowedParabolaNode;
            if (swallowedParabolaNode.Parent == null)
                return;

            Node leftAncestor, rightAncestor;
            swallowedParabolaNode.FindDirectionalAncestors(out leftAncestor, out rightAncestor);

            var leftAncestorEdgeRayData = (EdgeRayData) leftAncestor.Data;
            var rightAncestorEdgeRayData = (EdgeRayData) rightAncestor.Data;

            state.VoronoiEdges.Add(new VoronoiEdge(
                leftAncestorEdgeRayData.Origin, leftAncestorEdgeRayData.IsOriginBounded, e.Circumcenter, true));
            state.VoronoiEdges.Add(new VoronoiEdge(
                rightAncestorEdgeRayData.Origin, leftAncestorEdgeRayData.IsOriginBounded, e.Circumcenter, true));

            var leftParabolaNode = leftAncestor.Left.GetRightmost();
            var leftParabolaFocus = ((ParabolaData) leftParabolaNode.Data).Focus;

            var rightParabolaNode = rightAncestor.Right.GetLeftmost();
            var rightParabolaFocus = ((ParabolaData) rightParabolaNode.Data).Focus;

            state.DelanayEdges.Add(new VoronoiEdge(leftParabolaFocus, true, rightParabolaFocus, true));

            // One ancestor will become the edge shared by the left/right parabolas,
            // the other will be deleted. It's guranteed our parent is one of these
            // ancestors and that the other ancestor is above it (that is, it has a
            // parent) so opt to delete our parent.
            var nodeParent = swallowedParabolaNode.Parent;
            var nodeSibling = nodeParent.Left == swallowedParabolaNode ? nodeParent.Right : nodeParent.Left;
            NodeOperations.ReplaceNode(nodeParent, nodeSibling, ref state.Root);
            nodeParent.Left = nodeParent.Right = null;

            var olderAncestor = nodeParent == leftAncestor ? rightAncestor : leftAncestor;
            var leftFocusRightFocus = rightParabolaFocus - leftParabolaFocus; // x is positive
            var edgeDirection = new TVector2(-leftFocusRightFocus.Y, leftFocusRightFocus.X);
            var leftRightFocusCenter = new TVector2(
                MathUtil.Average(leftParabolaFocus.X, rightParabolaFocus.X),
                MathUtil.Average(leftParabolaFocus.Y, rightParabolaFocus.Y));
            var dy = sweepY - leftRightFocusCenter.Y;
            var edgeStart = e.Circumcenter;
            olderAncestor.Data = new EdgeRayData(true, edgeStart, edgeDirection);

            // add new potential circle events
            Node leftLeftLeaf;
            if (leftParabolaNode.TryGetLeftLeaf(out leftLeftLeaf))
                HandleAddCircleEvent(state, leftLeftLeaf, leftParabolaNode, rightParabolaNode, sweepY);

            Node rightRightLeaf;
            if (rightParabolaNode.TryGetRightLeaf(out rightRightLeaf))
                HandleAddCircleEvent(state, leftParabolaNode, rightParabolaNode, rightRightLeaf, sweepY);
        }
    }

    public class FortunesAlgorithmRenderer
    {
        private static readonly Pen _highlightPen = new Pen(Color.Red, 3);
        private static readonly Pen _passivePen = new Pen(Color.FromArgb(64, Color.Black), 1);

        // unhighlighted: active: peter river, else wet asphalt
        private static readonly Brush _unhighlightedInactiveParabolaBrush = new SolidBrush(Color.FromArgb(32, 52, 73, 94));
        private static readonly Brush _unhighlightedActiveParabolaBrush = new SolidBrush(Color.FromArgb(52, 152, 219));

        // highlighted: active: carrot, else emerald
        private static readonly Brush _highlightedInactiveParabolaBrush = new SolidBrush(Color.FromArgb(46, 204, 113));
        private static readonly Brush _highlightedActiveParabolaBrush = new SolidBrush(Color.Gold);

        private static readonly Pen _inProgressVoronoiEdgePen = new Pen(Color.FromArgb(52, 73, 94));
        private static readonly Pen _settledVoronoiEdgePen = new Pen(Color.FromArgb(44, 62, 80));
        private static readonly Pen _delaunayTriangulationEdgePen = new Pen(Color.FromArgb(231, 76, 60));

        private static readonly Pen _nextCircleEventPen = new Pen(Color.FromArgb(243, 156, 18), 3);
        private static readonly Pen _futureCircleEventPen = new Pen(Color.FromArgb(64, 243, 156, 18), 1);

        private readonly Font _font;

        public FortunesAlgorithmRenderer(Padding padding, Size boardSize, Font font)
        {
            Padding = padding;
            BoardSize = boardSize;
            _font = font;

            ImageSize = new Size(boardSize.Width + padding.Horizontal, boardSize.Height + padding.Vertical);
        }

        public Padding Padding { get; set; }
        public Size BoardSize { get; set; }
        public Size ImageSize { get; set; }

        public Bitmap Render(FortunesAlgorithmState state, TNumber? sweepY = null)
        {
            state.Root.ValidateTree();
            var bitmap = new Bitmap(ImageSize.Width, ImageSize.Height, PixelFormat.Format24bppRgb);

            using (var g = Graphics.FromImage(bitmap))
            {
                g.CompositingQuality = CompositingQuality.HighQuality;
                g.InterpolationMode = InterpolationMode.HighQualityBicubic;
                g.SmoothingMode = SmoothingMode.AntiAlias;

                g.Clear(Color.FromArgb(236, 240, 241));
                g.TranslateTransform(Padding.Left, Padding.Top);

                foreach (var initialPoint in state.InitialPoints)
                {
                    DrawPointSquare(g, initialPoint, 3);
                }

                foreach (var edge in state.VoronoiEdges)
                {
                    var start = edge.Start;
                    var end = edge.End;
                    var direction = end - start;

                    if (!edge.IsStartBounded)
                    {
                        start = end - direction * (end.Y + Padding.Top);
                    }

                    if (!edge.IsEndBounded)
                        end += direction * 100;

                    g.DrawLine(_settledVoronoiEdgePen, start.X, start.Y, end.X, end.Y);
                }

                foreach (var edge in state.DelanayEdges)
                {
                    var start = edge.Start;
                    var end = edge.End;
                    g.DrawLine(_delaunayTriangulationEdgePen, start.X, start.Y, end.X, end.Y);
                }

                Node highlightedNode = null;
                var eventQueueCopy = state.EventQueue.Copy();
                ISweepEvent sweepEvent;
                var eventCounter = 0;
                while (eventQueueCopy.TryDequeue(out sweepEvent))
                {
                    var eventPoint = sweepEvent.Point;
//                    g.DrawString(eventCounter.ToString(), _font, Brushes.Black, eventPoint.X, eventPoint.Y);

                    if (sweepEvent is CircleEvent)
                    {
                        var circleEvent = (CircleEvent) sweepEvent;

                        if (eventCounter == 0) highlightedNode = circleEvent.SwallowedParabolaNode;
                        var pen = eventCounter == 0 ? _nextCircleEventPen : _futureCircleEventPen;
                        DrawPointX(g, eventPoint, 3, pen);

                        g.DrawEllipse(
                            pen,
                            circleEvent.Circumcenter.X - circleEvent.Radius,
                            circleEvent.Circumcenter.Y - circleEvent.Radius,
                            circleEvent.Radius * 2,
                            circleEvent.Radius * 2
                        );
                    }
                    else
                    {
                        var siteEvent = (SiteEvent) sweepEvent;
                        if (state.Root != null && sweepY.HasValue && eventCounter == 0)
                        {
                            var parabolaNode = state.Root.FindDeepestNodeAtX(siteEvent.Point.X, siteEvent.Y);
                            var parabolaData = (ParabolaData) parabolaNode.Data;
                            var parabolaPoint = MathUtil.ComputeParabolaPointGivenX(
                                siteEvent.Point.X, sweepY.Value, parabolaData.Focus);
                            g.DrawLine(Pens.Black, parabolaPoint.X, parabolaPoint.Y, siteEvent.Point.X, siteEvent.Point.Y);
                        }
                    }
                    eventCounter++;
                }

                if (sweepY.HasValue)
                {
                    DrawBeachlineNode(g, state.Root, sweepY.Value, highlightedNode);
                    g.DrawLine(Pens.Black, -Padding.Left, (int)sweepY, -Padding.Left + ImageSize.Width, (int)sweepY);
                }
            }

            return bitmap;
        }

        private void DrawPointX(Graphics g, TVector2 point, int dims, Pen pen = null)
        {
            g.DrawLine(pen ?? Pens.Black, point.X - dims, point.Y - dims, point.X + dims, point.Y + dims);
            g.DrawLine(pen ?? Pens.Black, point.X + dims, point.Y - dims, point.X - dims, point.Y + dims);
        }

        private void DrawPointSquare(Graphics g, TVector2 point, int dims)
        {
            var offset = dims / 2;
            g.FillRectangle(
                Brushes.Black,
                point.X - offset,
                point.Y - offset,
                dims, dims
            );
        }

        private void DrawBeachlineNode(Graphics g, Node node, TNumber sweepY, Node highlightedNode = null)
        {
            if (node == null) return;

            if (node.IsLeaf())
            {
                var parabolaData = (ParabolaData) node.Data;

                TNumber renderStartX = -Padding.Left, renderEndX = BoardSize.Width + Padding.Right;
                TNumber activeStartX = 0, activeEndX = BoardSize.Width;
                Node leftEdgeRayNode, rightEdgeRayNode;
                node.FindDirectionalAncestors(out leftEdgeRayNode, out rightEdgeRayNode);
                if (leftEdgeRayNode != null)
                {
                    var leftEdgeRayData = (EdgeRayData) leftEdgeRayNode.Data;
                    var intersect = MathUtil.ComputeRayParabolaIntersection(
                        leftEdgeRayData.Origin, leftEdgeRayData.Direction,
                        parabolaData.Focus, sweepY);
                    activeStartX = intersect.X;
                    g.DrawLine(_inProgressVoronoiEdgePen, leftEdgeRayData.Origin.X, leftEdgeRayData.Origin.Y, intersect.X, intersect.Y);
                }
                if (rightEdgeRayNode != null)
                {
                    var rightEdgeRayData = (EdgeRayData) rightEdgeRayNode.Data;
                    var intersect = MathUtil.ComputeRayParabolaIntersection(
                        rightEdgeRayData.Origin, rightEdgeRayData.Direction,
                        parabolaData.Focus, sweepY);
                    activeEndX = intersect.X;
                    g.DrawLine(_inProgressVoronoiEdgePen, rightEdgeRayData.Origin.X, rightEdgeRayData.Origin.Y, intersect.X, intersect.Y);
                }
                DrawParabola(g, renderStartX, renderEndX, activeStartX, activeEndX, parabolaData.Focus, sweepY, node == highlightedNode);
            }
            else
            {
                DrawBeachlineNode(g, node.Left, sweepY, highlightedNode);
                DrawBeachlineNode(g, node.Right, sweepY, highlightedNode);
            }
        }

        private void DrawParabola(
            Graphics g, TNumber startX, TNumber endX, TNumber activeStartX, TNumber activeEndX,
            TVector2 focus, TNumber sweepY, bool isHighlighted)
        {
//                Console.WriteLine($"Drawing Parabola ({focus}, {sweepY}): [{startX}, {endX}) [{activeStartX}, {activeEndX})");

            DrawPointSquare(g, focus, 3);

            if (MathUtil.ApproximatelyEqual(focus.Y, sweepY))
                sweepY += 0.01f;

            for (int i = 0; i + startX < endX; i++)
            {
                TNumber x = startX + i;
                var p = MathUtil.ComputeParabolaPointGivenX(x, sweepY, focus);
                var point = new PointF((float)p.X, (float)p.Y);
                var activeBrush = isHighlighted ? _highlightedActiveParabolaBrush : _unhighlightedActiveParabolaBrush;
                var passiveBrush = isHighlighted ? _highlightedInactiveParabolaBrush : _unhighlightedInactiveParabolaBrush;
                if (activeStartX <= x && x <= activeEndX)
                    g.FillRectangle(activeBrush, point.X - 1, point.Y - 1, 3, 3);
                else
                    g.FillRectangle(passiveBrush, point.X, point.Y, isHighlighted ? 3 : 1, isHighlighted ? 3 : 1);
            }
        }

        public static FortunesAlgorithmRenderer Create(IReadOnlyCollection<Vector2> sites)
        {
            var minX = sites.Select(s => s.X).Min();
            var minY = sites.Select(s => s.Y).Min();
            var maxX = sites.Select(s => s.X).Max();
            var maxY = sites.Select(s => s.Y).Max();

            maxX += minX;
            maxY += minY;
            minX = 0;
            minY = 0;

            return new FortunesAlgorithmRenderer(
                new Padding(200, 100, 200, 100),
                new Size((int) Math.Ceiling(maxX - minX), (int) Math.Ceiling(maxY - minY)),
                SystemFonts.DefaultFont);
        }
    }

    /// <summary>
    /// Note: Order matters - circle should be ordered before site events.
    /// </summary>
    public enum SweepEventType
    {
        Circle = 0,
        Site = 1
    }

    public interface ISweepEvent {
        TVector2 Point { get; }
        TNumber Y { get; }
        SweepEventType Type { get; }
    }

    public class SiteEvent : ISweepEvent
    {
        public SiteEvent(TVector2 point)
        {
            Point = point;
        }

        public TVector2 Point { get; }
        public TNumber Y => Point.Y;
        public SweepEventType Type => SweepEventType.Site;

        public override string ToString() => $"{Type} {Point}";
    }

    public class CircleEvent : ISweepEvent
    {
        public CircleEvent(TVector2 circumcenter, TNumber radius, Node swallowedParabolaNode)
        {
            Circumcenter = circumcenter;
            Radius = radius;
            SwallowedParabolaNode = swallowedParabolaNode;
        }

        public TVector2 Circumcenter { get; set; }
        public TNumber Radius { get; set; }
        public TVector2 Point => Circumcenter + new TVector2(0, Radius);
        public TNumber Y => Point.Y;
        public Node SwallowedParabolaNode { get; }
        public SweepEventType Type => SweepEventType.Circle;

        public override string ToString() => $"{Type} {Y}; {((ParabolaData) SwallowedParabolaNode.Data).Focus}";
    }

    public class ImageDisplay : Form
    {
        private const int trackbarHeight = 50;

        private readonly List<Bitmap> frames = new List<Bitmap>();
        private readonly TrackBar trackbar;
        private readonly PictureBox display;

        public ImageDisplay(Size imageSize)
        {
            trackbar = new TrackBar
            {
                Orientation = Orientation.Horizontal,
                Location = new Point(0, imageSize.Height),
                Size = new Size(imageSize.Width, trackbarHeight),
                Minimum = 0,
                Maximum = 0,
                Value = 0,
                Enabled = false
            };
            display = new PictureBox
            {
                Location = new Point(0, 0),
                ClientSize = imageSize
            };
            Controls.Add(trackbar);
            Controls.Add(display);
            ClientSize = new Size(imageSize.Width, imageSize.Height + trackbarHeight);

            trackbar.ValueChanged += (s, e) => HandleDisplayFrame(trackbar.Value);
        }

        public void AddFrame(Bitmap frame)
        {
            frames.Add(frame);
            BeginInvoke(new Action(HandleFrameAdded));
        }

        private void HandleFrameAdded()
        {
            trackbar.Maximum = frames.Count - 1;
            trackbar.Enabled = true;
            HandleDisplayFrame(trackbar.Value);
        }

        private void HandleDisplayFrame(int index)
        {
            Text = $"Frame {index + 1} / {frames.Count}";
            display.Image = frames[index];
        }

        public static ImageDisplay CreateAndRunInBackground(Size imageSize)
        {
            ImageDisplay imageDisplay = null;
            var initializedEvent = new ManualResetEvent(false);
            var thread = new Thread(() =>
            {
                imageDisplay = new ImageDisplay(imageSize);
                imageDisplay.Shown += (s, e) => initializedEvent.Set();
                Application.Run(imageDisplay);
            });
            thread.SetApartmentState(ApartmentState.STA);
            thread.IsBackground = false;
            thread.Start();
            initializedEvent.WaitOne();
            return imageDisplay;
        }
    }

    public class VoronoiEdge
    {
        public VoronoiEdge(Vector2 start, bool isStartBounded, Vector2 end, bool isEndBounded)
        {
            Start = start;
            IsStartBounded = isStartBounded;
            End = end;
            IsEndBounded = isEndBounded;
        }

        public TVector2 Start { get; }
        public bool IsStartBounded { get; }

        public TVector2 End { get; }
        public bool IsEndBounded { get; }
    }

    public class EdgeRayData
    {
        public EdgeRayData(bool isOriginBounded, TVector2 origin, TVector2 direction)
        {
            IsOriginBounded = isOriginBounded;
            Origin = origin;
            Direction = direction;
        }

        public bool IsOriginBounded { get; }
        public TVector2 Origin { get; }
        public TVector2 Direction { get; }

        public override string ToString() => $"Edge {Origin}; {Direction}";
    }

    public class ParabolaData
    {
        public ParabolaData(TVector2 focus)
        {
            Focus = focus;
        }

        public TVector2 Focus { get; }
        public override string ToString() => $"Parabola {Focus}";
    }

    public class Program
    {
        public static void Main(string[] args)
        {
            var points = new HashSet<Vector2>
            {
                new Vector2(100, 100),
                new Vector2(200, 100),
                new Vector2(120, 150),
                new Vector2(180, 150)
            };

//            var scale = 1f;
//            var points = new HashSet<Vector2>
//            {
//                new Vector2(100 * scale, 100 * scale),
//                new Vector2(250 * scale, 100 * scale),
//                new Vector2(160 * scale, 150 * scale),
//                new Vector2(180 * scale, 160 * scale)
//            };

//            var random = new Random(1);
//            points = new HashSet<Vector2>(
//                Enumerable.Range(0, 30)
//                    .Select(i => new TVector2((TNumber) random.NextDouble() * 400, (TNumber) random.NextDouble() * 400))
//            );

            points = new HashSet<TVector2>(
                from i in Enumerable.Range(0, 4)
                from j in Enumerable.Range(0, 4)
                let shift = i % 2 == 0 ? 25 : 0
                let spacing = 50
                select new TVector2(j * spacing + shift, i * spacing)
            );

            var renderer = FortunesAlgorithmRenderer.Create(points);
            var display = ImageDisplay.CreateAndRunInBackground(renderer.ImageSize);

            var state = FortunesAlgorithmState.CreateAndInitialize(points);
            var fortunesAlgorithmExecutor = new FortunesAlgorithmExecutor();

//            // show end result
//            fortunesAlgorithmExecutor.ExecuteToCompletion(state);
//            display.AddFrame(renderer.Render(state));

            // animate
            var frames = new ObservableCollection<Bitmap>();
            frames.CollectionChanged += (s, e) => display.AddFrame((Bitmap)e.NewItems[0]);

            TNumber sweepY = -renderer.Padding.Top;
            ISweepEvent sweepEvent;
            while (state.EventQueue.TryPeek(out sweepEvent))
            {
                while (sweepY < sweepEvent.Y)
                {
                    TNumber stepSize = 1;
                    var bottomY = renderer.BoardSize.Height + renderer.Padding.Bottom;
                    if (sweepY > bottomY)
                    {
                        var speed = (TNumber)Math.Pow(sweepY - bottomY, 1.01);
                        stepSize = Math.Max(speed, stepSize);
                    }
                    frames.Add(renderer.Render(state, sweepY));
                    sweepY += stepSize;
                }

                fortunesAlgorithmExecutor.ExecuteNextEvent(state);
                frames.Add(renderer.Render(state, sweepEvent.Y + 0.1f));

                if (state.EventQueue.Count <= 5)
                {
                    fortunesAlgorithmExecutor.ExecuteToCompletion(state);
                }
            }
            while (sweepY < renderer.BoardSize.Height + renderer.Padding.Bottom)
                frames.Add(renderer.Render(state, sweepY++));

            fortunesAlgorithmExecutor.ProcessCleanup(state);
            frames.Add(renderer.Render(state));

            if (!Directory.Exists("output"))
                Directory.CreateDirectory("output");

            int frame = 0;
            for (int i = 0; i < 10; i++)
                frames.Last().Save($"output/{frame++}.gif", ImageFormat.Gif);
            for (int i = frames.Count - 1; i >= 0; i -= 14)
                frames[i].Save($"output/{frame++}.gif", ImageFormat.Gif);
            for (int i = 0; i < frames.Count; i += 3)
                frames[i].Save($"output/{frame++}.gif", ImageFormat.Gif);

//            using (var outputGifStream = File.OpenWrite("output.gif"))
//            using (var gifEncoder = new GifEncoder(
//                    outputGifStream, renderer.ImageSize.Width, renderer.ImageSize.Height))
//            {
//                gifEncoder.AddFrame(frames.Last(), TimeSpan.FromMilliseconds(1000));
//                for (int i = frames.Count - 1; i >= 0; i -= 10)
//                    gifEncoder.AddFrame(frames[i], TimeSpan.FromMilliseconds(10));
//                for (int i = 0; i < frames.Count; i++)
//                    gifEncoder.AddFrame(frames[i], TimeSpan.FromMilliseconds(10));
//            }
//            using (MagickImageCollection collection = new MagickImageCollection())
//            {
//                collection.Add(new MagickImage(frames.Last()));
//                collection[collection.Count - 1].AnimationDelay = 1000;
//
//                for (int i = frames.Count - 1; i >= 0; i -= 10)
//                {
//                    collection.Add(new MagickImage(frames[i]));
//                    collection[collection.Count - 1].AnimationDelay = 10;
//                }
//                for (int i = 0; i < frames.Count; i++)
//                {
//
//                    collection.Add(new MagickImage(frames[i]));
//                    collection[collection.Count - 1].AnimationDelay = 10;
//                }
//            }
        }
    }
}