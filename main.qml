
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

ApplicationWindow {
    visible: true
    width: 1280; height: 800
    color: "#05060A"

    // Background Gradient
    Rectangle {
        anchors.fill: parent
        gradient: Gradient {
            GradientStop { position: 0.0; color: "#1a1c2e" }
            GradientStop { position: 1.0; color: "#090A10" }
        }
    }

    RowLayout {
        anchors.fill: parent
        spacing: 0

        // Sidebar
        Rectangle {
            width: 200; height: parent.height
            color: "#151722"
            border.color: "#2a2d3e"
            
            Column {
                anchors.top: parent.top; anchors.topMargin: 40
                width: parent.width; spacing: 20
                
                Text { text: "APEX SYSTEM"; color: "#00E5FF"; font.bold: true; anchors.horizontalCenter: parent.horizontalCenter }
                
                Button {
                    text: "Import Model"; width: 160; anchors.horizontalCenter: parent.horizontalCenter
                    onClicked: backend.importStl("test.stl")
                }
            }
        }

        // Main Workspace
        Item {
            Layout.fillWidth: true; Layout.fillHeight: true
            
            // Glass Panel
            Rectangle {
                anchors.centerIn: parent
                width: 900; height: 600
                color: Qt.rgba(1, 1, 1, 0.05)
                radius: 20
                border.color: Qt.rgba(1, 1, 1, 0.1)
                
                Text {
                    anchors.centerIn: parent
                    text: backend.status
                    color: "#FFFFFF"
                    font.pixelSize: 32
                    font.weight: Font.Light
                }
            }
        }
    }
}
