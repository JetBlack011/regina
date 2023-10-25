//
//  ReginaDocument.swift
//  Regina
//
//  Created by Ben Burton on 23/10/2023.
//  Copyright © 2023 Regina Development Team. All rights reserved.
//

import SwiftUI
import UniformTypeIdentifiers

extension UTType {
    static var reginaData: UTType {
        UTType(importedAs: "org.computop.regina-data")
    }
}

struct ReginaDocument: FileDocument {
    var text: String

    init(text: String = "Hello, world!") {
        self.text = text
    }

    static var readableContentTypes: [UTType] { [.reginaData] }

    init(configuration: ReadConfiguration) throws {
        guard let data = configuration.file.regularFileContents,
              let string = String(data: data, encoding: .utf8)
        else {
            throw CocoaError(.fileReadCorruptFile)
        }
        text = string
    }
    
    func fileWrapper(configuration: WriteConfiguration) throws -> FileWrapper {
        let data = text.data(using: .utf8)!
        return .init(regularFileWithContents: data)
    }
}