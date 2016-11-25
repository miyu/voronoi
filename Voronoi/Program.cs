using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Drawing.Imaging;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices.ComTypes;
using System.Runtime.Remoting.Messaging;
using System.Threading;
using System.Windows.Forms;
using System.Windows.Forms.VisualStyles;

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
            const TFloat epsilon = 1E-5f;
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

    public class FortunesAlgorithmExecutor
    {
        private readonly PriorityQueue<ISweepEvent> _eventQueue = new PriorityQueue<ISweepEvent>(Comparer<ISweepEvent>.Create(CompareSweepEvents));
        private Node _root;

        public void Execute(HashSet<TVector2> inputs)
        {
            foreach (var input in inputs)
                _eventQueue.Enqueue(new SiteEvent(input));

            var display = Display.CreateAndRunInBackground(inputs);

            ISweepEvent sweepEvent;
            int lastFrame = 0;
            TNumber sweepY = 0;
            while (_eventQueue.TryPeek(out sweepEvent))
            {
//                while (sweepY < sweepEvent.Y)
//                {
//                    display.AddFrame(_root, _eventQueue, sweepY);
//                    sweepY++;
//                }
                sweepY = sweepEvent.Y;

                display.AddFrame(_root, _eventQueue, sweepEvent.Y);
                _eventQueue.TryDequeue(out sweepEvent);

                Console.WriteLine($"[{++lastFrame}] Processing Event: " + sweepEvent);
                if (sweepEvent.Type == SweepEventType.Site)
                    ProcessSiteEvent((SiteEvent) sweepEvent);
                else
                    ProcessCircleEvent((CircleEvent) sweepEvent);
            }
            for (int i = 1; i < 30; i++)
                display.AddFrame(_root, _eventQueue, sweepY + 0.1f * i);

//            while (sweepY < 800)
//            {
//                display.AddFrame(_root, _eventQueue, sweepY);
//                sweepY++;
//            }
            new CountdownEvent(1).Wait();
        }

        private void PrintBeachlineTree(Node highlightedNode = null)
        {
            PrintBeachlineTreeHelper(_root, highlightedNode, 1);
            Console.WriteLine("");
        }

        private void PrintBeachlineTreeHelper(Node current, Node highlightedNode, int indent)
        {
            if (current == null) return;
            Console.WriteLine(new string('\t', indent) + current.Data + (highlightedNode == current ? " ***" : ""));
            PrintBeachlineTreeHelper(current.Left, highlightedNode, indent + 1);
            PrintBeachlineTreeHelper(current.Right, highlightedNode, indent + 1);
        }

        private void ProcessSiteEvent(SiteEvent e)
        {
            PrintBeachlineTree();

            var newParabolaFocus = e.Point;
            if (_root == null)
            {
                _root = new Node(new ParabolaData(newParabolaFocus));
                PrintBeachlineTree(_root);
                return;
            }

            var cutNode = _root.FindDeepestNodeAtX(newParabolaFocus.X, newParabolaFocus.Y);
            if (!cutNode.IsLeaf()) throw new NotSupportedException("Cutting edge node not supported.");

            var cutParabolaFocus = ((ParabolaData)cutNode.Data).Focus;
            if (cutParabolaFocus.Y == newParabolaFocus.Y)
            {
                Node leftParabola, rightParabola;
                Node newNode = BeachlineNodeOperations.ComputeSameYParabolaNodeCut(
                    cutParabolaFocus, newParabolaFocus,
                    out leftParabola, out rightParabola);
                NodeOperations.ReplaceNode(cutNode, newNode, ref _root);

                Node leftLeftParabola;
                if (leftParabola.TryGetLeftLeaf(out leftLeftParabola))
                    HandleAddCircleEvent(leftLeftParabola, leftParabola, rightParabola);

                Node rightRightParabola;
                if (rightParabola.TryGetRightLeaf(out rightRightParabola))
                    HandleAddCircleEvent(leftParabola, rightParabola, rightRightParabola);

                PrintBeachlineTree(newNode);
            }
            else
            {
                Node leftParabola, centerParabola, rightParabola;
                Node newNode = BeachlineNodeOperations.ComputeThreeParabolasFromDifferentYParabolaNodeCut(
                    cutParabolaFocus, newParabolaFocus,
                    out leftParabola, out centerParabola, out rightParabola);
                NodeOperations.ReplaceNode(cutNode, newNode, ref _root);

                Node leftLeftParabola;
                if (leftParabola.TryGetLeftLeaf(out leftLeftParabola))
                    HandleAddCircleEvent(leftLeftParabola, leftParabola, centerParabola);

                Node rightRightParabola;
                if (rightParabola.TryGetRightLeaf(out rightRightParabola))
                    HandleAddCircleEvent(centerParabola, rightParabola, rightRightParabola);

                PrintBeachlineTree(newNode);
            }
        }

        void HandleAddCircleEvent(Node leftParabola, Node centerParabola, Node rightParabola)
        {
            var leftParabolaFocus = ((ParabolaData) leftParabola.Data).Focus;
            var centerParabolaFocus = ((ParabolaData) centerParabola.Data).Focus;
            var rightParabolaFocus = ((ParabolaData) rightParabola.Data).Focus;

            // sign of left.y - center.y
            var leftCenterSign = leftParabolaFocus.Y.CompareTo(centerParabolaFocus.Y);

            // sign of right.y - center.y
            var rightCenterSign = rightParabolaFocus.Y.CompareTo(centerParabolaFocus.Y);
            if ((leftCenterSign > 0 && rightCenterSign >= 0) ||
                (rightCenterSign > 0 && leftCenterSign >= 0))
            {
                // center gets swallowed by left/right when sweep line hits bottom of foci circumcircle
                TVector2 circumcenter;
                TNumber radius;
                MathUtil.FindCircumcircle(
                    leftParabolaFocus, centerParabolaFocus, rightParabolaFocus,
                    out circumcenter, out radius);
                _eventQueue.Enqueue(
                    new CircleEvent(
                        circumcenter,
                        radius,
                        centerParabola
                    ));
            }
        }

        void ProcessCircleEvent(CircleEvent e)
        {
            PrintBeachlineTree(e.SwallowedParabolaNode);

            var sweepY = e.Y;
            var swallowedParabolaNode = e.SwallowedParabolaNode;
            if (swallowedParabolaNode.Parent == null)
                return;

            Node leftAncestor, rightAncestor;
            swallowedParabolaNode.FindDirectionalAncestors(out leftAncestor, out rightAncestor);

            var leftParabolaNode = leftAncestor.Left.GetRightmost();
            var leftParabolaFocus = ((ParabolaData) leftParabolaNode.Data).Focus;

            var rightParabolaNode = rightAncestor.Right.GetLeftmost();
            var rightParabolaFocus = ((ParabolaData) rightParabolaNode.Data).Focus;

            // One ancestor will become the edge shared by the left/right parabolas,
            // the other will be deleted. It's guranteed our parent is one of these
            // ancestors and that the other ancestor is above it (that is, it has a
            // parent) so opt to delete our parent.
            var nodeParent = swallowedParabolaNode.Parent;
            var nodeSibling = nodeParent.Left == swallowedParabolaNode ? nodeParent.Right : nodeParent.Left;
            NodeOperations.ReplaceNode(nodeParent, nodeSibling, ref _root);
            nodeParent.Left = nodeParent.Right = null;

            var olderAncestor = nodeParent == leftAncestor ? rightAncestor : leftAncestor;
            var leftFocusRightFocus = rightParabolaFocus - leftParabolaFocus; // x is positive
            var edgeDirection = new TVector2(-leftFocusRightFocus.Y, leftFocusRightFocus.X);
            var leftRightFocusCenter = new TVector2(
                MathUtil.Average(leftParabolaFocus.X, rightParabolaFocus.X),
                MathUtil.Average(leftParabolaFocus.Y, rightParabolaFocus.Y));
            var dy = sweepY - leftRightFocusCenter.Y;
            var edgeStart = e.Circumcenter;
            olderAncestor.Data = new EdgeRayData(edgeStart, edgeDirection);

            // add new potential circle events
            Node leftLeftLeaf;
            if (leftParabolaNode.TryGetLeftLeaf(out leftLeftLeaf))
                HandleAddCircleEvent(leftLeftLeaf, leftParabolaNode, rightParabolaNode);

            Node rightRightLeaf;
            if (rightParabolaNode.TryGetRightLeaf(out rightRightLeaf))
                HandleAddCircleEvent(leftParabolaNode, rightParabolaNode, rightRightLeaf);

            PrintBeachlineTree(e.SwallowedParabolaNode);
        }

        public static int CompareSweepEvents(ISweepEvent a, ISweepEvent b)
        {
            var result = a.Y.CompareTo(b.Y);
            if (result != 0) return result;
            return ((int) a.Type).CompareTo((int) b.Type);
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

        public class Display : Form
        {
            private static readonly Pen highlightPen = new Pen(Color.Red, 3);
            private const int trackbarHeight = 50;

            private readonly Padding padding;
            private readonly Size boardSize;
            private readonly Size paddedImageSize;

            private readonly List<Bitmap> frames = new List<Bitmap>();
            private TrackBar trackbar;
            private PictureBox display;

            public Display(Padding padding, Size boardSize)
            {
                this.padding = padding;
                this.boardSize = boardSize;
                paddedImageSize = new Size(padding.Horizontal + boardSize.Width, padding.Vertical + boardSize.Height);

                InitializeComponents();
            }

            private void InitializeComponents()
            {
                trackbar = new TrackBar
                {
                    Orientation = Orientation.Horizontal,
                    Location = new Point(0, paddedImageSize.Height),
                    Size = new Size(paddedImageSize.Width, trackbarHeight),
                    Minimum = 0,
                    Maximum = 0,
                    Value = 0,
                    Enabled = false
                };
                display = new PictureBox
                {
                    Location = new Point(0, 0),
                    ClientSize = paddedImageSize
                };
                Controls.Add(trackbar);
                Controls.Add(display);
                ClientSize = new Size(paddedImageSize.Width, paddedImageSize.Height + trackbarHeight);

                trackbar.ValueChanged += (s, e) => HandleDisplayFrame(trackbar.Value);
            }

            public void AddFrame(Node beachlineRoot, PriorityQueue<ISweepEvent> eventQueue, TNumber? sweepY)
            {
                var frame = RenderAlgorithmState(beachlineRoot, eventQueue, sweepY);
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

            private Bitmap RenderAlgorithmState(Node beachlineRoot, PriorityQueue<ISweepEvent> eventQueue, TNumber? sweepY)
            {
                beachlineRoot.ValidateTree();
                var bitmap = new Bitmap(paddedImageSize.Width, paddedImageSize.Height, PixelFormat.Format32bppRgb);

                using (var g = Graphics.FromImage(bitmap))
                {
                    g.Clear(Color.Lime);
                    g.TranslateTransform(padding.Left, padding.Top);

                    g.DrawRectangle(Pens.Black, 0, 0, boardSize.Width, boardSize.Height);

                    var eventQueueCopy = eventQueue.Copy();
                    ISweepEvent sweepEvent;
                    var eventCounter = 0;
                    while (eventQueueCopy.TryDequeue(out sweepEvent))
                    {
                        var eventPoint = sweepEvent.Point;
                        var pen = eventCounter == 0 ? highlightPen : Pens.Black;
                        DrawPointX(g, eventPoint, 3, pen);
                        g.DrawString(eventCounter.ToString(), Font, Brushes.Black, eventPoint.X, eventPoint.Y);
                        if (sweepEvent is CircleEvent)
                        {
                            var circleEvent = (CircleEvent) sweepEvent;
                            g.DrawEllipse(
                                pen,
                                circleEvent.Circumcenter.X - circleEvent.Radius,
                                circleEvent.Circumcenter.Y - circleEvent.Radius,
                                circleEvent.Radius * 2,
                                circleEvent.Radius * 2
                            );
                        }
                        eventCounter++;
                    }

                    if (sweepY.HasValue)
                    {
                        DrawBeachlineNode(g, beachlineRoot, sweepY.Value);
                        g.DrawLine(Pens.Black, -padding.Left, (int)sweepY, -padding.Left + paddedImageSize.Width, (int)sweepY);
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

            private void DrawBeachlineNode(Graphics g, Node node, TNumber sweepY)
            {
                if (node == null) return;

                if (node.IsLeaf())
                {
                    var parabolaData = (ParabolaData) node.Data;

                    TNumber renderStartX = -padding.Left, renderEndX = boardSize.Width + padding.Right;
                    TNumber highlightStartX = 0, highlightEndX = boardSize.Width;
                    Node leftEdgeRayNode, rightEdgeRayNode;
                    node.FindDirectionalAncestors(out leftEdgeRayNode, out rightEdgeRayNode);
                    if (leftEdgeRayNode != null)
                    {
                        var leftEdgeRayData = (EdgeRayData) leftEdgeRayNode.Data;
                        var intersect = MathUtil.ComputeRayParabolaIntersection(
                            leftEdgeRayData.Origin, leftEdgeRayData.Direction,
                            parabolaData.Focus, sweepY);
                        highlightStartX = intersect.X;
                        g.DrawLine(Pens.Red, leftEdgeRayData.Origin.X, leftEdgeRayData.Origin.Y, intersect.X, intersect.Y);

//                        var leftParabola = leftEdgeRayNode.Left.GetRightmost();
//                        var leftParabolaData = (ParabolaData) leftParabola.Data;
//                        Vector2 leftIntersection, rightIntersection;
//                        MathUtil.ComputeParabolaIntersections(
//                            parabolaData.Focus, leftParabolaData.Focus,
//                            sweepY, out leftIntersection, out rightIntersection
//                        );
//                        var originToLeftIntersect = leftIntersection - leftEdgeRayData.Origin;
//                        var originToRightIntersect = rightIntersection - leftEdgeRayData.Origin;
//                        originToLeftIntersect /= originToLeftIntersect.Length();
//                        originToRightIntersect /= originToRightIntersect.Length();
//                        var rayDirectionUnit = leftEdgeRayData.Direction / leftEdgeRayData.Direction.Length();
//                        var leftCross = rayDirectionUnit.X * originToLeftIntersect.Y -
//                                        rayDirectionUnit.Y * originToLeftIntersect.X;
//                        var rightCross = rayDirectionUnit.X * originToRightIntersect.Y -
//                                         rayDirectionUnit.Y * originToRightIntersect.X;
//                        highlightStartX = Math.Abs(leftCross) < Math.Abs(rightCross)
//                            ? leftIntersection.X
//                            : rightIntersection.X;
                        Console.WriteLine("L " + leftEdgeRayData);
                    }
                    if (rightEdgeRayNode != null)
                    {
                        var rightEdgeRayData = (EdgeRayData) rightEdgeRayNode.Data;
                        var intersect = MathUtil.ComputeRayParabolaIntersection(
                            rightEdgeRayData.Origin, rightEdgeRayData.Direction,
                            parabolaData.Focus, sweepY);
                        highlightEndX = intersect.X;
                        g.DrawLine(Pens.Red, rightEdgeRayData.Origin.X, rightEdgeRayData.Origin.Y, intersect.X, intersect.Y);
                        Console.WriteLine("R " + rightEdgeRayData);
                    }
                    DrawParabola(g, renderStartX, renderEndX, highlightStartX, highlightEndX, parabolaData.Focus, sweepY);
                }
                else
                {
                    DrawBeachlineNode(g, node.Left, sweepY);
                    DrawBeachlineNode(g, node.Right, sweepY);
                }
            }

            private void DrawParabola(
                Graphics g, TNumber startX, TNumber endX, TNumber highlightStartX, TNumber highlightEndX,
                TVector2 focus, TNumber sweepY)
            {
                Console.WriteLine($"Drawing Parabola ({focus}, {sweepY}): [{startX}, {endX}) [{highlightStartX}, {highlightEndX})");

                DrawPointSquare(g, focus, 3);

                if (MathUtil.ApproximatelyEqual(focus.Y, sweepY))
                    sweepY += 0.01f;

                for (int i = 0; i + startX < endX; i++)
                {
                    TNumber x = startX + i;
                    var p = MathUtil.ComputeParabolaPointGivenX(x, sweepY, focus);
                    var point = new PointF((float)p.X, (float)p.Y);
                    if (highlightStartX <= x && x <= highlightEndX)
                        g.FillRectangle(Brushes.Red, point.X - 1, point.Y - 1, 3, 3);
                    else
                        g.FillRectangle(Brushes.Black, point.X, point.Y, 1, 1);
                }
            }

            public static Display CreateAndRunInBackground(IReadOnlyCollection<TVector2> sites)
            {
                var minX = sites.Select(s => s.X).Min();
                var minY = sites.Select(s => s.Y).Min();
                var maxX = sites.Select(s => s.X).Max();
                var maxY = sites.Select(s => s.Y).Max();

                maxX += minX;
                maxY += minY;
                minX = 0;
                minY = 0;

                return CreateAndRunInBackground(
                    new Padding(100),
                    new Size((int) Math.Ceiling(maxX - minX), (int) Math.Ceiling(maxY - minY)));
            }

            public static Display CreateAndRunInBackground(Padding padding, Size boardSize)
            {
                Display display = null;
                var initializedEvent = new ManualResetEvent(false);
                var thread = new Thread(() =>
                {
                    display = new Display(padding, boardSize);
                    display.Shown += (s, e) => initializedEvent.Set();
                    Application.Run(display);
                });
                thread.SetApartmentState(ApartmentState.STA);
                thread.IsBackground = false;
                thread.Start();
                initializedEvent.WaitOne();
                return display;
            }
        }
    }

    public class EdgeRayData
    {
        public EdgeRayData(TVector2 origin, TVector2 direction)
        {
            Origin = origin;
            Direction = direction;
        }

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
//            var points = new HashSet<Vector2>
//            {
//                new Vector2(100, 100),
//                new Vector2(200, 100),
//                new Vector2(120, 150),
//                new Vector2(180, 150)
//            };

            var scale = 1f;
            var points = new HashSet<Vector2>
            {
                new Vector2(100 * scale, 100 * scale),
                new Vector2(250 * scale, 100 * scale),
                new Vector2(160 * scale, 150 * scale),
                new Vector2(180 * scale, 160 * scale)
            };

//            var random = new Random(1);
//            points = new HashSet<Vector2>(
//                Enumerable.Range(0, 5)
//                    .Select(i => new TVector2((TNumber) random.NextDouble() * 400, (TNumber) random.NextDouble() * 400))
//            );

            var fortunesAlgorithmExecutor = new FortunesAlgorithmExecutor();
            fortunesAlgorithmExecutor.Execute(points);
        }
    }
}