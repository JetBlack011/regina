
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Swift User Interface                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2023, Ben Burton                                   *
 *  For further details contact Ben Burton (bab@debian.org).              *
 *                                                                        *
 *  This program is free software; you can redistribute it and/or         *
 *  modify it under the terms of the GNU General Public License as        *
 *  published by the Free Software Foundation; either version 2 of the    *
 *  License, or (at your option) any later version.                       *
 *                                                                        *
 *  As an exception, when this program is distributed through (i) the     *
 *  App Store by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or     *
 *  (iii) Google Play by Google Inc., then that store may impose any      *
 *  digital rights management, device limits and/or redistribution        *
 *  restrictions that are required by its terms of service.               *
 *                                                                        *
 *  This program is distributed in the hope that it will be useful, but   *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *  General Public License for more details.                              *
 *                                                                        *
 *  You should have received a copy of the GNU General Public             *
 *  License along with this program; if not, write to the Free            *
 *  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,       *
 *  MA 02110-1301, USA.                                                   *
 *                                                                        *
 **************************************************************************/

import SwiftUI
import ReginaEngine

// TODO: We need to BAN held() for shared packets. Deep copies are leading to dangling pointers/references.

// TODO: Work out what parts of this interface need to be made scrollable for large links.

extension regina.StrandRefAlt: Identifiable {
    public var id: Int { id() }
}

extension regina.SharedLink {
    func strandsForComponent(index: Int) -> [regina.StrandRefAlt] {
        let start = component(index)

        if start.isNull() {
            return []
        }
        
        var ans = [regina.StrandRefAlt]()
        var s = start
        repeat {
            ans.append(s)
            s = s.next()
        } while !(s == start)
        return ans
    }
}

/**
 * A class that supports refreshing a view even if the underlying link has not changed.
 *
 * An example of where this is useful is, for instance, when the Jones polynomial has been
 * computed (and so there is new information to display, even though the link itself did not change).
 */
class ObservedLink: ObservableObject {
    let packet: regina.SharedLink
    
    init(packet: regina.SharedLink) {
        self.packet = packet
    }
    
    func refresh() {
        objectWillChange.send()
    }
}

enum LinkTab {
    case Crossings, Polynomials, Algebra, Codes, Graphs
}

struct LinkView: View {
    let packet: regina.SharedLink
    @State private var selection: LinkTab = .Crossings
    @Environment(\.horizontalSizeClass) var sizeClass
    
    init(packet: regina.SharedLink) {
        self.packet = packet
    }

    var body: some View {
        let link = packet.held()
        
        VStack {
            switch link.size() {
            case 0:
                switch link.countComponents() {
                case 0:
                    Text("Empty link")
                case 1:
                    if sizeClass == .compact {
                        Text("Unknot")
                        Text("No crossings")
                    } else {
                        Text("Unknot with no crossings")
                    }
                default:
                    if sizeClass == .compact {
                        Text("Unlink")
                        Text("\(link.countComponents()) components")
                        Text("No crossings")
                    } else {
                        Text("Unlink with \(link.countComponents()) components, no crossings")
                    }
                }
            case 1:
                // This must be alternating, and must have ≥1 component.
                let signs = (link.writhe() > 0 ? Text("+")
                    .foregroundColor(Color("Positive")) : Text("−")
                    .foregroundColor(Color("Negative")))
                if link.countComponents() == 1 {
                    if sizeClass == .compact {
                        Text("Alternating knot")
                        Text("1 crossing (\(signs))")
                    } else {
                        Text("Alternating knot with 1 crossing (\(signs))")
                    }
                } else {
                    if sizeClass == .compact {
                        Text("Alternating link")
                        Text("\(link.countComponents()) components")
                        Text("1 crossing (\(signs))")
                    } else {
                        Text("Alternating link with \(link.countComponents()) components, 1 crossing (\(signs))")
                    }
                }
            default:
                let pos = (link.writhe() + link.size()) / 2
                let neg = link.size() - pos
                let signs = (neg == 0 ? Text("all +").foregroundColor(Color("Positive")) : pos == 0 ? Text("all −").foregroundColor(Color("Negative")) : Text("\(pos) +").foregroundColor(Color("Positive")) + Text(", ") + Text("\(neg) −").foregroundColor(Color("Negative")))
    
                let alt = (link.isAlternating() ? "Alternating" : "Non-alternating")
                if link.countComponents() == 1 {
                    if sizeClass == .compact {
                        Text("\(alt) knot")
                        Text("\(link.size()) crossings (\(signs))")
                    } else {
                        Text("\(alt) knot with \(link.size()) crossings (\(signs))")
                    }
                } else {
                    if sizeClass == .compact {
                        Text("\(alt) link")
                        Text("\(link.countComponents()) components")
                        Text("\(link.size()) crossings (\(signs))")
                    } else {
                        Text("\(alt) link with \(link.countComponents()) components, \(link.size()) crossings (\(signs))")
                    }
                }
            }
            
            // TODO: Ipad 11" portrait, tab icons jump around when selected??
            // TODO: Make a persistent default tab
            TabView(selection: $selection) {
                LinkCrossingsView(packet: packet).tabItem {
                    Image(selection == .Crossings ? "Tab-Crossings-Bold" : "Tab-Crossings").renderingMode(.template)
                    Text("Crossings")
                }.tag(LinkTab.Crossings)
                LinkPolynomialsView(packet: packet).tabItem {
                    Image("Tab-Polynomials").renderingMode(.template)
                    Text("Polynomials")
                }.tag(LinkTab.Polynomials)
                LinkAlgebraView(packet: packet).tabItem {
                    Image("Tab-Algebra").renderingMode(.template)
                    Text("Algebra")
                }.tag(LinkTab.Algebra)
                LinkCodesView(packet: packet).tabItem {
                    Image("Tab-Codes").renderingMode(.template)
                    Text("Codes")
                }.tag(LinkTab.Codes)
                LinkGraphsView(packet: packet).tabItem {
                    Image(selection == .Graphs ? "Tab-Graphs-Bold" : "Tab-Graphs").renderingMode(.template)
                    Text("Graphs")
                }.tag(LinkTab.Graphs)
            }
        }
    }
}

struct LinkCrossingsView: View {
    let packet: regina.SharedLink
    // TODO: Make a persistent default display type
    @State private var pictures: Bool = true

    // TODO: Check that these update when switching dark/light mode.
    static private var posColour = Color("Positive")
    static private var negColour = Color("Negative")

    @AppStorage("displayUnicode") private var unicode = true

    /**
     * Returns an image depicting the given crossing, without the crossing index.
     */
    func iconFor(_ s: regina.StrandRefAlt) -> Image {
        // TODO: Should we preload these?
        // TODO: Re-render these at a larger size. We use then at 34pt but they are currently only rendered at 22pt.
        if s.crossing().sign() > 0 {
            if s.strand() == 1 {
                return Image("Crossing+U")
            } else {
                return Image("Crossing+L")
            }
        } else {
            if s.strand() == 1 {
                return Image("Crossing-U")
            } else {
                return Image("Crossing-L")
            }
        }
    }
    
    func pictureFor(_ s: regina.StrandRefAlt) -> some View {
        let textHeight = fontSize(forTextStyle: .body)
        let iconSize = textHeight * 1.7
        return ZStack {
            iconFor(s).resizable().frame(width: iconSize, height: iconSize)
            // TODO: Work out the x offset properly, using the text width.
            Text("\(s.crossing().index())")
                .offset(x: textHeight, y: -textHeight)
        }
        // TODO: Choose the padding dynamically.
        .padding(.top, 10.0)
    }

    func textFor(_ s: regina.StrandRefAlt) -> Text {
        let c = s.crossing()
        let colour = (c.sign() > 0 ? LinkCrossingsView.posColour : LinkCrossingsView.negColour)
        
        if unicode {
            if c.sign() > 0 {
                if s.strand() == 1 {
                    return Text("\(c.index())⁺").foregroundColor(colour)
                } else {
                    return Text("\(c.index())₊").foregroundColor(colour)
                }
            } else {
                if s.strand() == 1 {
                    return Text("\(c.index())⁻").foregroundColor(colour)
                } else {
                    return Text("\(c.index())₋").foregroundColor(colour)
                }
            }
        } else {
            if c.sign() > 0 {
                if s.strand() == 1 {
                    return Text("\(c.index())^+").foregroundColor(colour)
                } else {
                    return Text("\(c.index())_+").foregroundColor(colour)
                }
            } else {
                if s.strand() == 1 {
                    return Text("\(c.index())^-").foregroundColor(colour)
                } else {
                    return Text("\(c.index())_-").foregroundColor(colour)
                }
            }
        }
    }
    
    var body: some View {
        VStack(alignment: .leading) {
            HStack {
                Spacer()
                Picker("Display crossings:", selection: $pictures) {
                    Text("Pictures").tag(true)
                    Text("Text").tag(false)
                }.fixedSize()
                Spacer()
            }
            .padding(.vertical)

            List {
                ForEach(0..<packet.countComponents(), id: \.self) { i in
                    Section("Component \(i)") {
                        let strands = packet.strandsForComponent(index: i)
                        if strands.isEmpty {
                            Text("Unknot, no crossings")
                        } else if pictures {
                            // TODO: Fix sizes: 45 is about right for a 17-point font size
                            LazyVGrid(columns: [.init(.adaptive(minimum: 45, maximum: 45))]) {
                                ForEach(strands) { s in
                                    pictureFor(s)
                                }
                            }
                        } else {
                            // TODO: Fix sizes: 25 is about right for a 17-point font size with single digits.
                            LazyVGrid(columns: [.init(.adaptive(minimum: 25, maximum: 25))]) {
                                ForEach(strands) { s in
                                    textFor(s)
                                }
                            }
                        }
                    }
                }
            }.listStyle(.plain)

            Spacer()
        }.padding(.horizontal).textSelection(.enabled)
    }
}

enum HomflyStyle {
    case az, lm
}

struct LinkPolynomialsView: View {
    static let maxAuto = 6;

    @ObservedObject var observed: ObservedLink
    // TODO: Make a persistent HOMFLY-PT selection
    @State private var homflyStyle: HomflyStyle = .az
    @AppStorage("displayUnicode") private var unicode = true

    init(packet: regina.SharedLink) {
        observed = ObservedLink(packet: packet)
    }
    
    var body: some View {
        let link = observed.packet.held()

        // TODO: Add a "copy plain text" option on long press
        // TODO: When computing, use progress trackers with cancellation
        // TODO: When computing, can we disable the buttons???
        VStack(alignment: .leading) {
            Text("Jones").font(.headline).padding(.vertical)
            if link.knowsJones() || link.size() <= LinkPolynomialsView.maxAuto {
                var jones = observed.packet.jones()
                // TODO: Make utf-8 configurable
                if jones.isZero() || jones.minExp() % 2 == 0 {
                    let _: Void = jones.scaleDown(2)
                    if unicode {
                        Text(swiftString(jones.utf8("𝑡")))
                    } else {
                        Text(swiftString(jones.str("t")))
                    }
                } else {
                    if unicode {
                        Text(swiftString(jones.utf8("√𝑡")))
                    } else {
                        Text(swiftString(jones.str("sqrt_t")))
                    }
                }
            } else {
                Button("Compute…", systemImage: "gearshape") {
                    observed.packet.jones()
                    observed.refresh()
                }
            }
            HStack {
                Text("HOMFLY-PT").font(.headline).padding(.vertical)
                Spacer()
                Picker("HOMFLY-PT type", selection: $homflyStyle) {
                    Text("(𝛼, 𝑧)").tag(HomflyStyle.az)
                    Text("(ℓ, 𝑚)").tag(HomflyStyle.lm)
                }.pickerStyle(.segmented).fixedSize().labelsHidden()
            }
            if link.knowsHomfly() || link.size() <= LinkPolynomialsView.maxAuto {
                if homflyStyle == .az {
                    if unicode {
                        Text(swiftString(observed.packet.homflyAZ().utf8("𝛼", "𝑧")))
                    } else {
                        Text(swiftString(observed.packet.homflyAZ().str("a", "z")))
                    }
                } else {
                    if unicode {
                        Text(swiftString(observed.packet.homflyLM().utf8("ℓ", "𝑚")))
                    } else {
                        Text(swiftString(observed.packet.homflyLM().str("l", "m")))
                    }
                }
            } else {
                Button("Compute…", systemImage: "gearshape") {
                    observed.packet.homflyAZ()
                    observed.refresh()
                }
            }
            Text("Kauffman bracket").font(.headline).padding(.vertical)
            if link.knowsBracket() || link.size() <= LinkPolynomialsView.maxAuto {
                if unicode {
                    Text(swiftString(observed.packet.bracket().utf8("𝐴")))
                } else {
                    Text(swiftString(observed.packet.bracket().str("A")))
                }
            } else {
                Button("Compute…", systemImage: "gearshape") {
                    observed.packet.bracket()
                    observed.refresh()
                }
            }
            Spacer()
        }.padding(.horizontal).textSelection(.enabled)
    }
}

struct LinkAlgebraView: View {
    static let maxSimp = 50
    static let maxRecognise = 50

    @ObservedObject var observed: ObservedLink
    @State var simplifiedGroup: regina.GroupPresentation?
    @AppStorage("displayUnicode") private var unicode = true

    init(packet: regina.SharedLink) {
        observed = ObservedLink(packet: packet)
    }

    var body: some View {
        let link = observed.packet.held()
        let autoSimp = (link.size() <= LinkAlgebraView.maxSimp)
        let group = simplifiedGroup ?? link.group(autoSimp)
        
        VStack(alignment: .leading) {
            HStack {
                Spacer()
                Text(link.countComponents() == 1 ? "Knot Group" : "Link Group").font(.headline).padding(.vertical)
                Spacer()
            }

            if !autoSimp {
                Text("Not automatically simplified").italic().padding(.bottom)
            }
            
            if group.countRelations() <= LinkAlgebraView.maxRecognise {
                let name = group.recogniseGroup(unicode)
                if name.length() > 0 {
                    Text("Name: \(swiftString(name))").padding(.bottom)
                }
            }

            // TODO: Do we want to headline the line headings here?
            let nGen = group.countGenerators()
            let alphabetic = (nGen <= 26)
            if nGen == 0 {
                Text("No generators").padding(.bottom)
            } else if nGen == 1 {
                Text("1 generator: a").padding(.bottom)
            } else if nGen == 2 {
                Text("2 generators: a, b").padding(.bottom)
            } else if alphabetic {
                let lastGen = String(UnicodeScalar(96 + Int(nGen))!)
                Text("\(nGen) generators: a … \(lastGen)").padding(.bottom)
            } else {
                Text("\(nGen) generators: g0 … g\(nGen - 1)").padding(.bottom)
            }
            
            let nRel = group.countRelations()
            if nRel == 0 {
                Text("No relations").padding(.bottom)
            } else {
                Text(nRel == 1 ? "1 relation:" : "\(nRel) relations:").padding(.bottom)
                // TODO: Should we put the relations inside a visible frame?
                // TODO: Should we be using a List or a ScrollView?
                List {
                    // We are using internal pointers within group, so ensure that
                    // group survives this entire block:
                    withExtendedLifetime(group) {
                        ForEach(0..<group.countRelations(), id: \.self) { i in
                            let rel = group.__relationUnsafe(i).pointee
                            if unicode {
                                Text(swiftString(rel.utf8(alphabetic)))
                            } else {
                                Text(swiftString(rel.str(alphabetic)))
                            }
                        }
                    }
                }
                .listStyle(.plain)
                //.padding(.horizontal)
                // TODO: Verify that the scrollable area works as it should
            }
            
            HStack {
                Spacer()
                Button("Try to simplify", systemImage: "rectangle.compress.vertical") {
                    var working = group
                    working.intelligentSimplify()
                    // TODO: If we could not simplify, inform the user and do not update
                    simplifiedGroup = working
                    // Currently Regina's links do not have a way to receive
                    // the simplified group, since link groups are not cached.
                }
                Spacer()
            }.padding(.vertical)
            
            Spacer()
        }.padding(.horizontal)
    }
}

enum LinkCode {
    case gauss, dt, signature, pd, jenkins
}

struct LinkCodesView: View {
    let packet: regina.SharedLink
    // TODO: Make a persistent default code
    @State private var selected: LinkCode = .gauss

    @ViewBuilder func onlyKnots(code: String, plural: Bool) -> some View {
        let capitalised = code.prefix(1).capitalized + code.dropFirst()
        let detail = "\(capitalised) \(plural ? "are" : "is") currently only available for knots, not multi-component links."
        
        if #available(macOS 14.0, iOS 17.0, *) {
            ContentUnavailableView {
                // Label("No \(code)", systemImage: "link")
                Label {
                    Text("No \(code)")
                } icon: {
                    // TODO: What size should we be using here?
                    // TODO: Render a large (64pt) version of this icon
                    Image("Link").renderingMode(.template).resizable().frame(width: 64, height: 64)
                }
            } description: {
                Text(detail)
            }
        } else {
            Text(detail)
        }
    }
    
    var body: some View {
        VStack(alignment: .leading) {
            HStack {
                Spacer()
                Picker("Display code:", selection: $selected) {
                    Text("Gauss codes").tag(LinkCode.gauss)
                    Text("Dowker-Thistlethwaite").tag(LinkCode.dt)
                    Text("Knot signature").tag(LinkCode.signature)
                    Text("Planar diagram code").tag(LinkCode.pd)
                    Text("Jenkins format").tag(LinkCode.jenkins)
                }.fixedSize()
                Spacer()
            }
            .padding(.vertical)
            
            let link = packet.held()
            
            switch (selected) {
            case .gauss:
                if link.countComponents() == 1 {
                    Text("Classical: ").font(.headline) + Text((swiftString(link.gauss())))
                    (Text("Oriented: ").font(.headline) + Text(swiftString(link.orientedGauss()))).padding(.top)
                } else {
                    Spacer()
                    onlyKnots(code: "Gauss codes", plural: true)
                    Spacer()
                }
            case .dt:
                if link.countComponents() == 1 {
                    if link.size() > 26 {
                        Text(swiftString(link.dt(false)))
                    } else {
                        Text(swiftString(link.dt(true)))
                        Text(swiftString(link.dt(false))).padding(.top)
                    }
                } else {
                    Spacer()
                    onlyKnots(code: "Dowker-Thistlethwaite notation", plural: false)
                    Spacer()
                }
            case .signature:
                if link.countComponents() == 1 {
                    Text(swiftString(link.knotSig(true, true)))
                } else {
                    Spacer()
                    onlyKnots(code: "knot signatures", plural: true)
                    Spacer()
                }
           case .pd:
                Text(swiftString(link.pd()))
                if link.countTrivialComponents() > 0 {
                    Text("This link includes one or more zero-crossing unknot components. These components are omitted entirely from the planar diagram code.").padding(.top)
                }
                if link.pdAmbiguous() {
                    Text("This link includes at least one component that consists entirely of over-crossings. The planar diagram code does not carry enough information to reconstruct the orientation of such a component.").padding(.top)
                }
            case .jenkins:
                Text(swiftString(link.jenkins()))
            }
            Spacer()
        }.padding(.horizontal).textSelection(.enabled)
    }
}

enum LinkGraph {
    case tree, nice
}

struct LinkGraphsView: View {
    let packet: regina.SharedLink
    // TODO: Make a persistent default graph type
    @State private var selected: LinkGraph = .tree

    // TODO: Ensure the graphs are visible in dark mode also.
    
    var body: some View {
        VStack {
            HStack {
                Spacer()
                Picker("Display graph:", selection: $selected) {
                    Text("Tree decomposition").tag(LinkGraph.tree)
                    Text("Nice tree decomposition").tag(LinkGraph.nice)
                }.fixedSize()
                Spacer()
            }
            .padding(.vertical)
            
            var tree = regina.TreeDecomposition(packet.held(), .Upper)
            if selected == .nice {
                let _ = tree.makeNice(nil)
            }
            if tree.size() == 1 {
                Text("1 bag, width \(tree.width())").padding(.bottom)
            } else {
                Text("\(tree.size()) bags, width \(tree.width())").padding(.bottom)
            }

            SvgView(cxxString: regina.svgUsingDot(tree.dot()))
            Spacer()
        }.padding(.horizontal).textSelection(.enabled)
   }
}

struct LinkView_Previews: PreviewProvider {
    static var previews: some View {
        LinkView(packet: regina.SharedLink(regina.ExampleLink.whitehead()))
    }
}