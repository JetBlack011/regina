
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  iOS User Interface                                                    *
 *                                                                        *
 *  Copyright (c) 1999-2013, Ben Burton                                   *
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

#import "NewDim2TriangulationController.h"
#import "PacketTreeController.h"
#import "dim2/dim2exampletriangulation.h"
#import "dim2/dim2triangulation.h"

@interface NewDim2TriangulationController ()
@property (weak, nonatomic) IBOutlet UISegmentedControl *types;
@property (weak, nonatomic) IBOutlet UIView *container;
@property (weak, nonatomic) NewPacketPageViewController *pages;
@end

@implementation NewDim2TriangulationController

- (void)viewDidLoad
{
    self.pages = static_cast<NewPacketPageViewController*>(self.childViewControllers.lastObject);
    [self.pages fillWithPages:@[@"newDim2TriEmpty", @"newDim2TriExample", @"newDim2TriSurface", @"newDim2TriIsosig"]
                 pageSelector:self.types
                   defaultKey:@"NewDim2TriangulationPage"];
}

- (IBAction)create:(id)sender
{
    regina::NPacket* ans = [self.pages create];
    if (ans) {
        self.spec.parent->insertChildLast(ans);
        [self.spec.tree viewPacket:ans];
        [self dismissViewControllerAnimated:YES completion:nil];
    }
}

- (IBAction)cancel:(id)sender
{
    [self dismissViewControllerAnimated:YES completion:nil];
}

@end

#pragma mark - Empty page

@implementation NewDim2TriangulationEmptyPage

- (regina::NPacket*)create
{
    regina::NPacket* ans = new regina::Dim2Triangulation();
    ans->setPacketLabel("2-D Triangulation");
    return ans;
}

@end

#pragma mark - Example triangulation

typedef regina::Dim2Triangulation* (*Dim2TriangulationCreator)();

/**
 * Represents a single option in the examples picker.
 */
@interface Dim2ExampleTriangulation : NSObject

@property (strong, nonatomic) NSString* name;
@property (assign, nonatomic) Dim2TriangulationCreator creator;

+ (id)exampleWithName:(NSString*)name creator:(Dim2TriangulationCreator)creator;
- (regina::Dim2Triangulation*)create;

@end

@implementation Dim2ExampleTriangulation

+ (id)exampleWithName:(NSString *)name creator:(Dim2TriangulationCreator)creator
{
    Dim2ExampleTriangulation* e = [[Dim2ExampleTriangulation alloc] init];
    if (e) {
        e.name = name;
        e.creator = creator;
    }
    return e;
}

- (regina::Dim2Triangulation *)create
{
    regina::Dim2Triangulation* ans = (*self.creator)();
    ans->setPacketLabel(self.name.UTF8String);
    return ans;
}

@end

#pragma mark - Example page

@interface NewDim2TriangulationExamplePage () <UIPickerViewDataSource, UIPickerViewDelegate> {
    NSArray* options;
}
@property (weak, nonatomic) IBOutlet UIPickerView *example;
@end

#define KEY_LAST_EXAMPLE @"NewDim2TriangulationExample"

@implementation NewDim2TriangulationExamplePage

- (void)viewDidLoad
{
    options = @[[Dim2ExampleTriangulation exampleWithName:@"Sphere (2 triangles)" creator:&regina::Dim2ExampleTriangulation::sphere],
                [Dim2ExampleTriangulation exampleWithName:@"Sphere (tetrahedron boundary)" creator:&regina::Dim2ExampleTriangulation::sphereTetrahedron],
                [Dim2ExampleTriangulation exampleWithName:@"Sphere (octahedron boundary)" creator:&regina::Dim2ExampleTriangulation::sphereOctahedron],
                [Dim2ExampleTriangulation exampleWithName:@"Disc" creator:&regina::Dim2ExampleTriangulation::disc],
                [Dim2ExampleTriangulation exampleWithName:@"Annulus" creator:&regina::Dim2ExampleTriangulation::annulus],
                [Dim2ExampleTriangulation exampleWithName:@"Möbius band" creator:&regina::Dim2ExampleTriangulation::mobius],
                [Dim2ExampleTriangulation exampleWithName:@"Torus" creator:&regina::Dim2ExampleTriangulation::torus],
                [Dim2ExampleTriangulation exampleWithName:@"Projective plane" creator:&regina::Dim2ExampleTriangulation::rp2],
                [Dim2ExampleTriangulation exampleWithName:@"Klein bottle" creator:&regina::Dim2ExampleTriangulation::kb]];
   
    self.example.dataSource = self;
    self.example.delegate = self;
    
    [self.example selectRow:[[NSUserDefaults standardUserDefaults] integerForKey:KEY_LAST_EXAMPLE] inComponent:0 animated:NO];
}

- (NSInteger)numberOfComponentsInPickerView:(UIPickerView *)pickerView
{
    return 1;
}

- (NSInteger)pickerView:(UIPickerView *)pickerView numberOfRowsInComponent:(NSInteger)component
{
    return options.count;
}

- (NSString *)pickerView:(UIPickerView *)pickerView titleForRow:(NSInteger)row forComponent:(NSInteger)component
{
    return [options[row] name];
}

- (void)pickerView:(UIPickerView *)pickerView didSelectRow:(NSInteger)row inComponent:(NSInteger)component
{
    [[NSUserDefaults standardUserDefaults] setInteger:[self.example selectedRowInComponent:0] forKey:KEY_LAST_EXAMPLE];
}

- (regina::NPacket *)create
{
    return [options[[self.example selectedRowInComponent:0]] create];
}

@end

#pragma mark - Surface page

@interface NewDim2TriangulationSurfacePage ()
@property (weak, nonatomic) IBOutlet UISegmentedControl *orbl;
@property (weak, nonatomic) IBOutlet UILabel *genusExpln;
@property (weak, nonatomic) IBOutlet UITextField *genus;
@property (weak, nonatomic) IBOutlet UITextField *punctures;
@end

@implementation NewDim2TriangulationSurfacePage

// TODO

- (regina::NPacket*)create
{
    regina::NPacket* ans = new regina::Dim2Triangulation();
    ans->setPacketLabel("2-D Triangulation");
    return ans;
}

@end

#pragma mark - Isosig page

@interface NewDim2TriangulationIsosigPage ()
@property (weak, nonatomic) IBOutlet UITextField *isosig;
@end

@implementation NewDim2TriangulationIsosigPage

- (IBAction)editingEnded:(id)sender {
    NewDim2TriangulationController* c = static_cast<NewDim2TriangulationController*>(self.parentViewController.parentViewController);
    [c create:sender];
}

- (regina::NPacket *)create
{
    std::string sig = [self.isosig.text stringByTrimmingCharactersInSet:[NSCharacterSet whitespaceCharacterSet]].UTF8String;
    if (sig.empty()) {
        UIAlertView* alert = [[UIAlertView alloc]
                              initWithTitle:@"Empty isomorphism signature"
                              message:@"Please type an isomorphism signature into the box provided."
                              delegate:nil
                              cancelButtonTitle:@"Close"
                              otherButtonTitles:nil];
        [alert show];
        return 0;
    }
    
    regina::Dim2Triangulation* t = regina::Dim2Triangulation::fromIsoSig(sig);
    if (! t) {
        UIAlertView* alert = [[UIAlertView alloc]
                              initWithTitle:@"Invalid isomorphism signature"
                              message:@"I could not interpret the given isomorphism signature."
                              delegate:nil
                              cancelButtonTitle:@"Close"
                              otherButtonTitles:nil];
        [alert show];
        return 0;
    }
    
    t->setPacketLabel(sig);
    return t;
}

@end
